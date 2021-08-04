mrstModule add ad-core ad-blackoil deckformat diagnostics mrst-gui ad-props incomp optimization network-models

 
%% Run EGG field simulation

% Run a realization as a reference model 
realization = 9; 
Settup_Egg_simulation 
 
wellSols_ref =  wellSols;
model_ref    = model;
states_ref   = states;
schedule_ref = schedule;
W_ref        = schedule_ref.control.W;

%% Creating the network

              
% Defining a well structure with only once well cell per well, to impose a
% Network that has only once well-cell per well.
W_ref_V2 = W_ref;
for i = 1:numel(W_ref_V2) % TODO maybe this part should be done more systematically.
    W_ref_V2(i).cells = W_ref_V2(i).cells([1]);
end

network_type = 'fd_postprocessor';
%network_type = 'fd_preprocessor';
%network_type  = 'injectors_to_producers';
%network_type = 'all_to_all';
% network_type = 'user_defined_edges';
% edges        = [1 9;2 9;2 10;3 9;3 11;...
%                 4 9; 4 10; 4 11; 4 12;...
%                 5 10; 5 12; 6 11; 7 11; 7 12; 8 12];
switch network_type
    case 'all_to_all'
        ntwkr =  Network(W_ref_V2,model_ref.G,...
                                'type',network_type);
    case 'injectors_to_producers'
        ntwkr =  Network(W_ref_V2,model_ref.G,...
                                 'type',network_type,...
                                 'injectors',[1:8],...
                                 'producers',[9:12]);
    case 'user_defined_edges'
        ntwkr =  Network(W_ref_V2,model_ref.G,...
                                 'type',network_type,...
                                 'edges',edges);
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
                                 'state_number',40);
    otherwise
        error('\nType of network: %s is not implemented\n', network_type);     
        
end



                     
% Ploting network
if any(strcmp(network_type,{'fd_preprocessor','fd_postprocessor'}))
% Ploting network with flow-diagnostics parameters

    TT = ntwkr.network.Edges.Transmissibility;
    pv = ntwkr.network.Edges.PoreVolume;
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

% We first create a simple MRST model
L= nthroot(sum(model_ref.operators.pv./model_ref.rock.poro)*25,3);
G = cartGrid([10, 1, numedges(ntwkr.network)], [L, L/5 ,L/5]*meter^3);
G = computeGeometry(G);

fluid =  model_ref.fluid;
rock = makeRock(G, 200*milli*darcy, 0.1);

gravity off
model = GenericBlackOilModel(G, rock, fluid);
model.gas=false;
model.OutputStateFunctions = {};


% Then we maps the Network into the MRST model
obj = NetworkModel(model,10, ntwkr.network,W_ref);
model = obj.model;
W     = obj.W;
indexs.faces = obj.Graph.Edges.Face_Indices;
indexs.cells = obj.Graph.Edges.Cell_Indices;


% Reorganizing the parameters indices 
% for lumping in th ModelParameter class
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
state0 = initState(G, W , 400*barsa,[0.2, 0.8]); 
dt = schedule.step.val;
schedule = simpleSchedule(dt(1:40), 'W', W);
schedule_0=schedule;

%% Define the parameter for the calibration
 prob = struct('model', model, 'schedule', schedule, 'state0', state0);                               
 parameters =  {};                          
 parameters{1} = ModelParameter(prob, 'name', 'conntrans','scaling','log','relativeLimits',[.01 20]);
 parameters{2} = ModelParameter(prob, 'name', 'porevolume','lumping',cell_lumping,'subset', cell_sub,'relativeLimits', [.01 20]);
 parameters{3} = ModelParameter(prob, 'name', 'transmissibility','lumping',faces_lumping,'subset', face_sub,'relativeLimits', [.01 20]);


 %% Simulating the initial DD model

weighting =  {'WaterRateWeight',  (150/day)^-1, ...
              'OilRateWeight',    (150/day)^-1, ...
              'BHPWeight',        (40*barsa)^-1};    
 
values = applyFunction(@(p)p.getParameterValue(prob), parameters);
% Initialize well productivity index, pore volume and transmisibility
values{1} =  10*values{1};
if any(strcmp(network_type,{'fd_preprocessor','fd_postprocessor'}))
     values{2} =  pv/10;
     values{3} =  TT;
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
  
  %% Optimization/Calibration/History matching
  
objh = @(p)evaluateMatch(p,obj,prob,parameters,states_ref);

[v, p_opt, history] = unitBoxBFGS(p0_fd, objh,'objChangeTol',  1e-8, 'maxIt', 25, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Simulating all simulation time
schedule = simpleSchedule(dt, 'W', W);
prob.schedule = schedule;
 
 [~,~,wellSols_opt] = evaluateMatch(p_opt,obj,prob,parameters, states_ref,'Gradient','none');
 [~,~,wellSols_0] = evaluateMatch(p0_fd,obj,prob,parameters, states_ref,'Gradient','none');

plotWellSols({wellSols_ref,wellSols_0,wellSols_opt},{schedule_ref.step.val,schedule.step.val,schedule.step.val})
legend('reference model','initial DD model','optimize DD model')

