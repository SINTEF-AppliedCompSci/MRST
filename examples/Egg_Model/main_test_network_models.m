mrstModule add ad-core ad-blackoil deckformat diagnostics mrst-gui ad-props incomp optimization network-models

 
%% Run EGG field simulation
Settup_Egg_simulation 
 
wellSols_ref =  wellSols;
model_ref    = model;
states_ref   = states;
schedule_ref = schedule;
W_ref        = schedule_ref.control.W;

%% Creating the network


               
% Defining a well structure with only once well cell per well, to impose a
% Network that has only once well cell per well.
W_ref_V2 = W_ref;
for i = 1:numel(W_ref_V2)
    W_ref_V2(i).cells = W_ref_V2(i).cells([1]);
end

%network_type = 'fd_postprocessor';
%network_type = 'fd_preprocessor';
%network_type  = 'injectors_to_producers';
network_type = 'all_to_all';

% Network derive by flow diagnostics
switch network_type
    case 'all_to_all'
        ntwkr =  Network(W_ref_V2,model_ref.G,...
                                'type',network_type);
    case 'injectors_to_producers'
        ntwkr =  Network(W_ref_V2,model_ref.G,...
                                 'injectors',[1:8],...
                                 'producers',[9:12],...
                                 'type',network_type);
    case 'fd_preprocessor'
        ntwkr =  Network(W_ref_V2,model_ref.G,...
                                 'problem',problem,...
                                 'type',network_type,...
                                 'flow_filter',1*stb/day);
    case 'fd_postprocessor'    
        ntwkr =  Network(W_ref_V2,model_ref.G,...
                                 'problem',problem,...
                                 'type',network_type,...
                                 'flow_filter',1*stb/day,...
                                 'state_number',40);
    otherwise
        error('\nType of network: %s is not implemented\n', network_type);     
        
end



                     
% Ploting flow-diagnostics parameters
if any(strcmp(network_type,{'fd_preprocessor','fd_postprocessor'}))

    TT = ntwkr.network.Edges.Transmissibility;
    pv = ntwkr.network.Edges.PoreVolume;

    figure, subplot(1,2,1)
    ntwkr.plotNetwork('NetworkLineWidth',10*TT/max(TT))
    title('Transmissibility')
    axis off

    subplot(1,2,2)
    ntwkr.plotNetwork('NetworkLineWidth',10*pv/max(pv))
    title('PoreVolume')
    axis off
else
    ntwkr.plotNetwork()
end

%% Creatting data driven model (TODO put this section inside Network Model)

L = 435;
G = cartGrid([10, 1, numedges(ntwkr.network)], [L, L/5 ,L/5]*meter^3);
G = computeGeometry(G);


fluid =  model_ref.fluid;
rock = makeRock(G, 200*milli*darcy, 0.1);

gravity off
model = GenericBlackOilModel(G, rock, fluid);
model.gas=false;
model.OutputStateFunctions = {};



%%
%[model,W,indexs] = createDDmodel_1(model,10,DD.Graph,W_ref);

obj = NetworkModel(model,10, ntwkr.network,W_ref);

model = obj.model;
W     = obj.W;
indexs.faces = obj.Graph.Edges.Face_Indices;
indexs.cells = obj.Graph.Edges.Cell_Indices;




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



%% Preparing parameters and scaling values for each one

if any(strcmp(network_type,{'fd_preprocessor','fd_postprocessor'}))
    Tr_boxlimits = [0.1*TT , 1.4*TT];
    Pv_boxlimits = [0.01*pv/10 , 2*pv/10];
else
    Tr_boxlimits = [];
    Pv_boxlimits = [];
end

               


%% Prepare the model for simulation.

model = model.validateModel();
state0 = initState(G, W , 400*barsa,[0.2, 0.8]); 
dt = schedule.step.val;
schedule = simpleSchedule(dt(1:40), 'W', W);
schedule_0=schedule;


 prob = struct('model', model, 'schedule', schedule, 'state0', state0);                               
 parameters =  {};                          
 parameters{1} = ModelParameter(prob, 'name', 'conntrans','relativeLimits',[.01 20]);
 parameters{2} = ModelParameter(prob, 'name', 'porevolume','lumping',cell_lumping,'subset', cell_sub,'relativeLimits', [.01 20]);
 parameters{3} = ModelParameter(prob, 'name', 'transmissibility','lumping',faces_lumping,'subset', face_sub,'relativeLimits', [.01 20]);


 %% Simulating the initial DD model

weighting =  {'WaterRateWeight',  (150/day)^-1, ...
              'OilRateWeight',    (150/day)^-1, ...
              'BHPWeight',        (400*barsa)^-1};    
 
values = applyFunction(@(p)p.getParameterValue(prob), parameters);
% scale values
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
  
  %% Optimization
  
objh = @(p)evaluateMatch(p,obj,prob,parameters,states_ref);

[v, p_opt, history] = unitBoxBFGS(p0_fd, objh,'objChangeTol',  1e-8, 'maxIt', 25, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Simulating all simulation time
schedule = simpleSchedule(dt, 'W', W);
prob.schedule = schedule;
 
 [~,~,wellSols_opt] = evaluateMatch(p_opt,obj,prob,parameters, states_ref,'Gradient','none');
 [~,~,wellSols_0] = evaluateMatch(p0_fd,obj,prob,parameters, states_ref,'Gradient','none');


plotWellSols({wellSols_ref,wellSols_0,wellSols_opt},{schedule_ref.step.val,schedule.step.val,schedule.step.val})
legend('reference model','initial DD model','optimize DD model')

PlotEggmodel_results_2