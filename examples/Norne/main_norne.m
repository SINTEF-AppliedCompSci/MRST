mrstModule add ad-core ad-blackoil deckformat diagnostics...
               mrst-gui ad-props incomp optimization...
               network-models example-suite linearsolvers 

%% Setting up the reference model           
exampleNorne_FD

wellSols_ref =  wellSols;
model_ref    = example.model;
states_ref   = states;
state0_ref   = example.state0;
schedule_ref = example.schedule;
W_ref        = schedule_ref.control.W;


%  %%  compute diagnostics and create a data driven model  DD

%% Creating the network

% Defining a well structure with only once well cell per well, to impose a
% Network that has only once well-cell per well.

W_ref_V2 = W_ref;
for i = 1:numel(W_ref_V2) % TODO maybe this part should be done more systematically.
    num_cells = numel(W_ref_V2(i).cells)
    %W_ref_V2(i).cells = W_ref_V2(i).cells([1 round(num_cells/2) num_cells]);
    W_ref_V2(i).cells = W_ref_V2(i).cells([round(num_cells/2)]);
end

network_type = 'fd_postprocessor';
%network_type = 'fd_preprocessor';
%network_type  = 'injectors_to_producers';
%network_type = 'all_to_all';
%network_type = 'user_defined_edges';

switch network_type
    case 'all_to_all'
        ntwkr =  Network(W_ref_V2,model_ref.G,...
                                'type',network_type);
    case 'injectors_to_producers'
        ntwkr =  Network(W_ref_V2,model_ref.G,...
                                 'type',network_type,...
                                 'injectors',[1:6],...
                                 'producers',[7:11]);
    case 'fd_preprocessor'
        % Network derive by flow diagnostics      
        ntwkr =  Network(W_ref_V2,model_ref.G,...
                                 'type',network_type,...
                                 'problem',problem,...
                                 'flow_filter',1*stb/day);
    case 'fd_postprocessor'
        % Network derive by flow diagnostics
        ntwkr =  Network(W_ref_V2,model_ref.G,...
                                 'type',network_type,...
                                 'problem',problem,...
                                 'flow_filter',1*stb/day,...
                                 'state_number',6);
    otherwise
        error('\nType of network: %s is not implemented\n', network_type);     
        
end



                     
% Ploting network
if any(strcmp(network_type,{'fd_preprocessor','fd_postprocessor'}))

    TT = ntwkr.network.Edges.Transmissibility;
    pv = ntwkr.network.Edges.PoreVolume;
% Ploting network with flow-diagnostics parameters
    figure, subplot(1,2,1);
    ntwkr.plotNetwork('NetworkLineWidth',10*TT/max(TT));
    title('Transmissibility');
    axis off

    subplot(1,2,2);
    ntwkr.plotNetwork('NetworkLineWidth',10*pv/max(pv));
    title('PoreVolume');
    axis off;
else
    ntwkr.plotNetwork()
    axis off;
end

%% Creatting data driven model

L= nthroot(sum(model_ref.operators.pv./model_ref.rock.poro)*25,3)  ;                            
G = cartGrid([10, 1, numedges(ntwkr.network)], [L, L/5 ,L/5]*meter^3);
G = computeGeometry(G);


fluid = example.model.fluid;
scaling = {'SWL', .1, 'SWCR', .2, 'SWU', .9, 'SOWCR', .1, 'KRW', .9, 'KRO', .8};

rock = makeRock(G, 1000*milli*darcy, 0.2);

gravity off
model = GenericBlackOilModel(G, rock, fluid,'gas', false);
model = imposeRelpermScaling(model, scaling{:});

model.OutputStateFunctions = {};


nc = 10;

NetModel = NetworkModel(model,nc,...
                   ntwkr.network,...
                   W_ref);

model = NetModel.model;
W     = NetModel.W;
indexs.faces = NetModel.Graph.Edges.Face_Indices;
indexs.cells = NetModel.Graph.Edges.Cell_Indices;




faces_lumping = zeros(size(model.operators.T));
cell_lumping  = zeros(size(model.operators.pv));
for i = 1:numel(indexs.faces)
    faces_lumping(indexs.faces{i}) = i;
    cell_lumping(indexs.cells{i})  = i;
end
 faces_lumping(find(faces_lumping==0))=NaN;
 cell_lumping(find(cell_lumping==0))=NaN;


cell_sub = find(~isnan(cell_lumping));
cell_lumping = cell_lumping(cell_sub);
face_sub = find(~isnan(faces_lumping));
faces_lumping = faces_lumping(face_sub);


%% Preparing parameters



% Prepare the model for simulation.

model = model.validateModel();
model.toleranceCNV = 1e-6;

nc =  model.G.cells.num;
nf =  numel(model.operators.T);
state0 = initState(model.G, W , 100*barsa,[0,1]); 
dt = schedule_ref.step.val;

 schedule_training = simpleSchedule(dt(1:6), 'W', W);
 
% select configuration for sensitivity computations
 prob_training  = struct('model', model, 'schedule', schedule_training , 'state0', state0);                               
config = {...%name      include     scaling    boxlims     lumping     subset       relativeLimits
          'porevolume',       1,   'linear',       [],  cell_lumping,   cell_sub,   [.001 4]
          'conntrans',        1,   'log',          [],          [],       [] ,      [.001 100]
          'transmissibility', 1,   'log'   ,       [],  faces_lumping, face_sub ,   [.001 100]
          'swl',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'swcr',             1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'swu',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'sowcr',            1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'krw',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'kro',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'sw',               1,   'linear',   [0 .6],  cell_lumping,   cell_sub,  [0 10]
          'pressure'          1,   'linear',       [],  cell_lumping,   cell_sub,  [.1 4]};
parameters = [];
for k = 1:size(config,1)
    if config{k, 2} ~= 0
        parameters = addParameter(parameters, prob_training, 'name',    config{k,1}, ...
                                                     'scaling', config{k,3}, ...
                                                     'boxLims', config{k,4}, ...
                                                     'lumping', config{k,5}, ...
                                                     'subset',  config{k,6},...
                                                     'relativeLimits',config{k,7});
    end
end
 

figure, plotCellData(model.G, state0.s(:,2)), colorbar, view(0,0);axis equal tight;  daspect([1,0.1,0.1])


 %% Simulation of the reference model and Data-driven model

         weighting =  {'WaterRateWeight',  (500/day)^-1, ...
                       'OilRateWeight',    (500/day)^-1, ...
                       'BHPWeight',        (50*barsa())^-1};
 
values = applyFunction(@(p)p.getParameterValue(prob_training), parameters);
% scale values
      values{2} =  0.5*values{2};
 if any(strcmp(network_type,{'fd_preprocessor','fd_postprocessor'}))
      values{1} =  pv/10;
      values{3} =  TT;
 end
u = cell(size(values));
for k = 1:numel(u)    
    u{k} = parameters{k}.scale(values{k});
end
p0_fd = vertcat(u{:});  
 
  obj = @(model, states, schedule, states_ref, tt, tstep, state) matchObservedOW(model, states, schedule, states_ref,...
            'computePartials', tt, 'tstep', tstep, weighting{:},'state',state,'from_states',false);
       
 [misfitVal_0,~,wellSols_0,states_0] = evaluateMatch(p0_fd,obj,prob_training,parameters, states_ref,'Gradient','none');          
 
 plotWellSols({wellSols_ref,wellSols_0},{schedule_ref.step.val,prob_training.schedule.step.val})
  



%% Optimization
  
objh = @(p)evaluateMatch(p,obj,prob_training,parameters,states_ref);

[v, p_opt, history] = unitBoxBFGS(p0_fd, objh,'objChangeTol',  1e-8, 'maxIt',30, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Simulating all simulation time
prob = prob_training;
prob.schedule = simpleSchedule(dt, 'W', W);
 
 [~,~,wellSols_opt] = evaluateMatch(p_opt,obj,prob,parameters, states_ref,'Gradient','none');
 [~,~,wellSols_0] = evaluateMatch(p0_fd,obj,prob,parameters, states_ref,'Gradient','none');


fh = plotWellSols({wellSols_ref,wellSols_0,wellSols_opt},{schedule_ref.step.val,prob.schedule.step.val,prob.schedule.step.val})
set(fh, 'name','Norne')
legend('reference model','initial DD model','optimize DD model')
                                       