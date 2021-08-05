mrstModule add network-models ad-core ad-props ad-blackoil optimization diagnostics

%ad-core ad-blackoil deckformat diagnostics mrst-gui ad-props incomp optimization ddmodel


%% Run EGG field simulation
%Setting_up_Olympus 
 
%% Defining reference observations

wellSols_ref = wellSols;
model_ref    = problem.SimulatorSetup.model;
states_ref   = states;
schedule_ref = problem.SimulatorSetup.schedule;
W_ref        = schedule_ref.control.W;
state0_ref   = problem.SimulatorSetup.state0;

%% Creating the network

% Defining a well structure with only once well cell per well, to impose a
% Network that has only once well-cell per well.

W_ref_V2 = W_ref;
for i = 1:numel(W_ref_V2) % TODO maybe this part should be done more systematically.
    W_ref_V2(i).cells = W_ref_V2(i).cells([7]);
end

%network_type = 'fd_postprocessor';
%network_type = 'fd_preprocessor';
network_type  = 'injectors_to_producers';
%network_type = 'all_to_all';
%network_type = 'user_defined_edges';

switch network_type
    case 'all_to_all'
        ntwkr =  Network(W_ref_V2,model_ref.G,...
                                'type',network_type);
    case 'injectors_to_producers'
        ntwkr =  Network(W_ref_V2,model_ref.G,...
                                 'type',network_type,...
                                 'injectors',[1:7],...
                                 'producers',[8:18]);
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
                                 'state_number',20);
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


% fluid =  model_ref.fluid;
% % fluid =  initSimpleADIFluid('phases', 'WO',...
% %                             'mu' , [ 0.398, 3.5]*centi*poise,...
% %                             'rho', [1020, 850]*kilogram/meter^3, ...
% %                             'n' , [ 2, 2]);
fluid = initSimpleADIFluid('phases', 'WO',... 
                           'mu' , [.3, 3]*centi*poise,...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n', [2 2]);
fluid .krPts  = struct('w', [0 0 1 1], 'ow', [0 0 1 1]);
scaling = {'SWL', .1, 'SWCR', .2, 'SWU', .9, 'SOWCR', .1, 'KRW', .9, 'KRO', .8};
c = 5e-5/barsa;
p_ref = 200*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);


rock = makeRock(G, 200*milli*darcy, 0.2);

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


%% Prepare the model for simulation.

model = model.validateModel();
state0 = initState(G, W , 200*barsa,[0.3, 0.7]); 
dt = schedule_ref.step.val;
schedule = simpleSchedule(dt(1:20), 'W', W);
schedule_0=schedule;


 prob = struct('model', model, 'schedule', schedule, 'state0', state0);                               
%  parameters =  {};                          
%  parameters{1} = ModelParameter(prob, 'name', 'conntrans','scaling','log','relativeLimits',[.01 20]);
%  parameters{2} = ModelParameter(prob, 'name', 'porevolume','lumping',cell_lumping,'subset', cell_sub,'scaling','log','relativeLimits', [.001 2]);
%  parameters{3} = ModelParameter(prob, 'name', 'transmissibility','lumping',faces_lumping,'subset', face_sub,'scaling','log','relativeLimits', [.001 100]);

 
nc =  model.G.cells.num;
nf =  numel(model.operators.T);
% select configuration for sensitivity computations
config = {...%name      include     scaling    boxlims     lumping     subset       relativeLimits
          'porevolume',       1,   'linear',       [],  cell_lumping,   cell_sub,   [.001 2]
          'conntrans',        1,   'log',          [],          [],       [] ,      [.001 20]
          'transmissibility', 1,   'log'   ,       [],  faces_lumping, face_sub ,   [.001 20]
          'swl',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'swcr',             1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'swu',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'sowcr',            1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'krw',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'kro',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'sw',               1,   'linear',   [0 .3],          [],     [],   [0.001 10]
          'pressure'          1,   'linear',       [],          [],     [],   [0.001 10]};
parameters = [];
for k = 1:size(config,1)
    if config{k, 2} ~= 0
        parameters = addParameter(parameters, prob, 'name',    config{k,1}, ...
                                                     'scaling', config{k,3}, ...
                                                     'boxLims', config{k,4}, ...
                                                     'lumping', config{k,5}, ...
                                                     'subset',  config{k,6},...
                                                     'relativeLimits',config{k,7});
    end
end

 %% Simulating the initial DD model

weighting =  {'WaterRateWeight',  (200/day)^-1, ...
              'OilRateWeight',    (20/day)^-1};    
 
values = applyFunction(@(p)p.getParameterValue(prob), parameters);
% scale values
%      values{1} =  10*values{1};
% if any(strcmp(network_type,{'fd_preprocessor','fd_postprocessor'}))
%      values{2} =  pv/10;
%      values{3} =  TT;
% end
u = cell(size(values));
for k = 1:numel(u)    
    u{k} = parameters{k}.scale(values{k});
end
p0_fd = vertcat(u{:});  
 
  obj = @(model, states, schedule, states_ref, tt, tstep, state) matchObservedOW(model, states, schedule, states_ref,...
            'computePartials', tt, 'tstep', tstep, weighting{:},'state',state,'from_states',false);
       
 [misfitVal_0,~,wellSols_0,states_0] = evaluateMatch(p0_fd,obj,prob,parameters, states_ref,'Gradient','none');          
 
 plotWellSols({wellSols_ref,wellSols_0},{schedule_ref.step.val,schedule_0.step.val})
  
 
%% Optimization
  
objh = @(p)evaluateMatch(p,obj,prob,parameters,states_ref);

[v, p_opt, history] = unitBoxBFGS(p0_fd, objh,'objChangeTol',  1e-8, 'maxIt', 50, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Simulating all simulation time
schedule = simpleSchedule(dt, 'W', W);
prob.schedule = schedule;
 
 [~,~,wellSols_opt] = evaluateMatch(p_opt,obj,prob,parameters, states_ref,'Gradient','none');
 [~,~,wellSols_0] = evaluateMatch(p0_fd,obj,prob,parameters, states_ref,'Gradient','none');


plotWellSols({wellSols_ref,wellSols_0,wellSols_opt},{schedule_ref.step.val,schedule.step.val,schedule.step.val})
legend('reference model','initial DD model','optimize DD model')
                                       