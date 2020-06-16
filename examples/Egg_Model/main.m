mrstModule add ad-core ad-blackoil deckformat diagnostics mrst-gui ad-props incomp optimization

 
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

%% Preparing parameters and scaling values for each one

WellIP = [];
cell = [];
well_index = []; 
levels = 1;        
for  i = 1:12
    for l = 1 : levels
        well_index = [well_index; i,l];
    end
    WellIP = [WellIP; W(i).WI(1)];
end


Tr_boxlimits = [0.1*TT ; 1.4*TT]';
Pv_boxlimits = [0.01*pv/10 ; 2*pv/10]';
well_boxlimits = [ 0.01*WellIP , ...
                  7*WellIP];
               

               
          well_IP = struct('name','conntrans',...
                                   'type','value',...
                                   'boxLims', well_boxlimits,...
                                   'distribution','general',...
                                   'Indx',well_index);


          transmisibility_conection = struct('name','transmissibility',...
                                   'type','value',...
                                   'boxLims', Tr_boxlimits ,...
                                   'distribution','connection',...
                                   'Indx',{indexs.faces});
                               

          porevolume_conection = struct('name','porevolume',...
                                   'type','value',...
                                   'boxLims',Pv_boxlimits,...
                                   'distribution','connection',...
                                   'Indx',{indexs.cells});                                                                
                               
parameters =  {};                          
parameters{1} = transmisibility_conection;
parameters{2}  = porevolume_conection ;
parameters{3}  = well_IP ;

%% Prepare the model for simulation.

model = model.validateModel();
state0 = initState(G, W , 400*barsa,[0.2, 0.8]); 
dt = schedule.step.val;


 %% Simulating the initial DD model

 schedule = simpleSchedule(dt(1:40), 'W', W);

         weighting =  {'WaterRateWeight',  1e5, ...
                       'OilRateWeight',    1e5, ...
                       'BHPWeight',        1e-5};
 schedule_0=schedule;
 
                  
  val{1} = TT;
  val{2} = pv/10;
  val{3} = WellIP;
  p0_fd = value2control(val,parameters);
 
  
 [misfitVal_0,gradient,wellSols_0,states_0] = Simulate_BFGS(p0_fd,parameters,model,schedule_0,state0, wellSols_ref,weighting,1);          
  plotWellSols({wellSols_ref,wellSols_0},{schedule_ref.step.val,schedule_0.step.val})
  
  %% Optimization
  
obj_scaling     = abs(misfitVal_0);      % objective scaling  
objh = @(p)Simulate_BFGS(p,parameters,model,schedule_0,state0,  wellSols_ref,weighting,obj_scaling);

[v, p_opt, history] = unitBoxBFGS(p0_fd, objh,'gradTol',             1e-2, ...
                                              'objChangeTol',        0.5e-3);


%% Simulating all simulation time
 schedule = simpleSchedule(dt, 'W', W);
 
 [misfitVal_opt,gradient_opt,wellSols_opt] = Simulate_BFGS(p_opt,parameters,model,schedule,state0, wellSols_ref,weighting,obj_scaling);

  [misfitVal_0,gradient_0,wellSols_0] = Simulate_BFGS(p0_fd,parameters,model,schedule,state0, wellSols_ref,weighting,obj_scaling);



plotWellSols({wellSols_ref,wellSols_0,wellSols_opt},{schedule_ref.step.val,schedule.step.val,schedule.step.val})
legend('reference model','initial DD model','optimize DD model')

PlotEggmodel_results_2