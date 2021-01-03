mrstModule add network-models ad-core ad-props ad-blackoil

%ad-core ad-blackoil deckformat diagnostics mrst-gui ad-props incomp optimization ddmodel


%% Run EGG field simulation
Setting_up_Olympus 
 
%% Defining reference observations

wellSols_ref =  wellSols;
model_ref    = model;
states_ref   = states;
schedule_ref = schedule;
W_ref        = schedule_ref.control.W;
state0_ref   = state0;

%%  Compute diagnostics 

 DD = WellPairNetwork(model_ref,schedule_ref,states_ref,state0_ref,wellSols_ref);
 DD =  DD.filter_wps(1*stb/day);
 DD.plotWellPairConnections()
 DD.plotWellPairsData('subplot',[6,8])
 
 
 % Initializing parameters
 for i =  1:numel(DD.wps)
    pv(i) = DD.wps{i}.volume;
    TT(i) = DD.wps{i}.Tr;
 end

%% Creatting data driven model

L = 435;
G = cartGrid([10, 1,numedges(DD.Graph)], [L, L/5 ,L/5]*meter^3);
G = computeGeometry(G);


%fluid =  model_ref.fluid;
fluid =  initSimpleADIFluid('phases', 'WO',...
                            'mu' , [ 0.398, 3.5]*centi*poise,...
                            'rho', [1020, 850]*kilogram/meter^3, ...
                            'n' , [ 2, 2]);
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
state0 = initState(G, W , 200*barsa,[0.2, 0.8]); 
dt = schedule.step.val;

 %% Simulating the initial DD model

 schedule = simpleSchedule(dt(1:20), 'W', W);

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
  