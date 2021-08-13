mrstModule add network-models ad-core ad-props ad-blackoil optimization diagnostics


%% Setting up Olympus field simulation
setting_up_olympus_model

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
    num_cells = numel(W_ref_V2(i).cells)
    W_ref_V2(i).cells = W_ref_V2(i).cells([1 round(num_cells/2) num_cells]);
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
                                 'state_number',1);
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

fluid = initSimpleADIFluid('phases', 'WO',... 
                           'mu' , [.3, 3]*centi*poise,...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n', [2 2]);
fluid .krPts  = struct('w', [0 0 1 1], 'ow', [0 0 1 1]);
scaling = {'SWL', .1, 'SWCR', .2, 'SWU', .9, 'SOWCR', .1, 'KRW', .9, 'KRO', .8};

rock = makeRock(G,500*milli*darcy, 0.2);

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
indices.faces = NetModel.Graph.Edges.Face_Indices;
indices.cells = NetModel.Graph.Edges.Cell_Indices;

[cell_lumping,cell_sub] = reorganizeIndices(indices.cells);
[faces_lumping,face_sub] = reorganizeIndices(indices.faces);

%% Prepare the model for simulation.

model = model.validateModel();
model.toleranceCNV = 1e-6;
state0 = initState(G, W , 215*barsa,[0, 1]); 
dt = schedule_ref.step.val;
schedule = simpleSchedule(dt(1:15), 'W', W);
schedule_0=schedule;
prob = struct('model', model, 'schedule', schedule, 'state0', state0);                               



%% Prepare the parameters to be calibrated

nc =  model.G.cells.num;
nf =  numel(model.operators.T);
% select configuration for sensitivity computations
config = {...%name      include     scaling    boxlims     lumping     subset       relativeLimits
          'porevolume',       1,   'linear',       [],  cell_lumping,   cell_sub,   [.001 2]
          'conntrans',        1,   'log',          [],          [],       [] ,      [.001 10]
          'transmissibility', 1,   'log'   ,       [],  faces_lumping, face_sub ,   [.001 4]
          'swl',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'swcr',             1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'swu',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'sowcr',            1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'krw',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'kro',              1,   'linear',       [],  ones(nc,1),       [],       [.5 2]
          'sw',               1,   'linear',   [0 .99], cell_lumping,   cell_sub,  [0 9]
          'pressure'          0,   'linear',       [],          [],     [],   [0.001 10]};
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

weighting =  {'WaterRateWeight',  (100/day)^-1, ...
              'OilRateWeight',    (10/day)^-1};    
 
values = applyFunction(@(p)p.getParameterValue(prob), parameters);
% scale values
      values{2} =  values{2};
if any(strcmp(network_type,{'fd_preprocessor','fd_postprocessor'}))
     values{1} =  pv/10;
     values{3} =  TT;
else
  rng(12345)
  values{10} =  rand(size(values{10}));  
end

u = cell(size(values));
for k = 1:numel(u)    
    u{k} = parameters{k}.scale(values{k});
end
p0_fd = vertcat(u{:});  
 
  obj = @(model, states, schedule, states_ref, tt, tstep, state) matchObservedOW(model, states, schedule, states_ref,...
            'computePartials', tt, 'tstep', tstep, weighting{:},'state',state,'from_states',false);
       
 [misfitVal_0,~,wellSols_0,states_0] = evaluateMatch(p0_fd,obj,prob,parameters, states_ref,'Gradient','none');          
 
 plotWellSols({wellSols_ref,wellSols_0},{schedule_ref.step.val,schedule_0.step.val})
  
 
%% Gradient based Optimization BFGS
  
objh = @(p)evaluateMatch(p,obj,prob,parameters,states_ref);

[v, p_opt, history] = unitBoxBFGS(p0_fd, objh,'objChangeTol',  1e-8, 'maxIt', 40, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Simulating based cased and calibrated model for comparison
schedule = simpleSchedule(dt, 'W', W);
prob.schedule = schedule;
 
 [~,~,wellSols_opt] = evaluateMatch(p_opt,obj,prob,parameters, states_ref,'Gradient','none');
 [~,~,wellSols_0] = evaluateMatch(p0_fd,obj,prob,parameters, states_ref,'Gradient','none');

fh = plotWellSols({wellSols_ref,wellSols_0,wellSols_opt},{schedule_ref.step.val,schedule.step.val,schedule.step.val});
set(fh, 'name','Olympus Network models, HM adjoin')
legend('reference model','initial network model','calibrated network model')
                                       