mrstModule add ad-core ad-blackoil deckformat diagnostics mrst-gui ad-props incomp optimization network-models

 
%% Run EGG field simulation
Settup_Egg_simulation 
 
wellSols_ref =  wellSols;
model_ref    = model;
states_ref   = states;
schedule_ref = schedule;
W_ref        = schedule_ref.control.W;

%%  Compute diagnostics 

 DD = WellPairNetwork(model_ref,schedule_ref,states_ref,state,wellSols_ref);
 DD =  DD.filter_wps(1*stb/day);
 DD.plotWellPairConnections()
 DD.plotWellPairsData('subplot',[4,4])
 
 
 
 % Initializing parameters
 for i =  1:numel(DD.wps)
    pv(i) = DD.wps{i}.volume;
    TT(i) = DD.wps{i}.Tr;
end

%% Creatting data driven model

L = 435;
G = cartGrid([10, 1, 16], [L, L/5 ,L/5]*meter^3);
G = computeGeometry(G);


fluid =  model_ref.fluid;
rock = makeRock(G, 1000*milli*darcy, 0.1);

gravity off
model = GenericBlackOilModel(G, rock, fluid);
model.gas=false;
model.OutputStateFunctions = {};

 
[model,W,indexs] = createDDmodel_1(model,10,DD.Graph,W_ref);

model.toleranceCNV = 1e-6;



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

WellIP = [];
well_index = []; 
levels = 1;   
well_index_2 = {};
for  i = 1:numel(W)
    for l = 1 : levels
        well_index = [well_index; i,l];
    end
    WellIP = [WellIP; W(i).WI(1)];
end


Tr_boxlimits = [0.1*TT ; 1.4*TT]';
Pv_boxlimits = [0.01*pv/10 ; 2*pv/10]';
well_boxlimits = [ 0.01*WellIP , ...
                  7*WellIP];
               


%% Prepare the model for simulation.

model = model.validateModel();
state0 = initState(G, W , 400*barsa,[0.2, 0.8]); 
dt = schedule.step.val;
schedule = simpleSchedule(dt(1:40), 'W', W);
schedule_0=schedule;



 prob = struct('model', model, 'schedule', schedule, 'state0', state0);                               
 parameters =  {};                          
 parameters{1} = ModelParameter(prob, 'name', 'conntrans','relativeLimits',[.01 10]);
 parameters{2} = ModelParameter(prob, 'name', 'porevolume','lumping',cell_lumping,'subset', cell_sub,'boxLims', Pv_boxlimits);
 parameters{3} = ModelParameter(prob, 'name', 'transmissibility','lumping',faces_lumping,'subset', face_sub,'boxLims', Tr_boxlimits );


 %% Simulating the initial DD model

weighting =  {'WaterRateWeight',  (150/day)^-1, ...
              'OilRateWeight',    (150/day)^-1, ...
              'BHPWeight',        (400*barsa)^-1};    
 
values = applyFunction(@(p)p.getParameterValue(prob), parameters);
% scale values
 values{1} =  WellIP;
 values{2} =  pv'/10;
 values{3} =  TT';
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

[v, p_opt, history] = unitBoxBFGS(p0_fd, objh,'objChangeTol',  1e-4, 'maxIt', 25, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Simulating all simulation time
 porb.schedule = simpleSchedule(dt, 'W', W);
 [misfitVal_opt,~,wellSols_opt] = evaluateMatch(p_opt,obj,prob,parameters, states_ref,'Gradient','none');

 [misfitVal_0,~,wellSols_0] = evaluateMatch(p0_fd,obj,prob,parameters, states_ref,'Gradient','none');


plotWellSols({wellSols_ref,wellSols_0,wellSols_opt},{schedule_ref.step.val,schedule.step.val,schedule.step.val})
legend('reference model','initial DD model','optimize DD model')

PlotEggmodel_results_2