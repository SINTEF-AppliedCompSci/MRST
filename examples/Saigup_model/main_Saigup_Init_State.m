mrstModule add ad-core ad-blackoil deckformat diagnostics mrst-gui ad-props incomp optimization kpn-ds example-suite linearsolvers 

 mrstVerbose('off')
% 
 load 'Saigup.mat'

% 
wellSols_ref =  wellSols;
model_ref    = example.model;
states_ref   = states;
state0_ref   = example.state0;
schedule_ref = problem.SimulatorSetup.schedule;
W_ref        = schedule_ref.control.W;


%  %%  compute diagnostics and create a data driven model  DD

%%  Compute diagnostics 

 DD = WellPairNetwork(model_ref,schedule_ref,states_ref,state0_ref,wellSols_ref);
 DD =  DD.filter_wps(150*stb/day);
 DD.plotWellPairConnections()

%%

%load 'Data_and_Diagnostic_SAIGUP.mat'


tt =  cumsum(schedule_ref.step.val);
for i = 1 : numel(DD.wps)
      subplot(4,5,i)
    plot(tt/day, -DD.wps{i}.data.BFPD*day/stb,'k')
     title(DD.wps{i}.wellpair_name);
end
figure
for i = 1 : numel(DD.wps)
    subplot(4,5,i)
    plot(tt/day, DD.wps{i}.data.s_avg,'k')
     title(DD.wps{i}.wellpair_name);
end


%% Two parameter per well pair
k =1;
for i =  1:numel(DD.wps)
    pv([k k+1]) = DD.wps{i}.volume;
    TT([k k+1]) = DD.wps{i}.Tr;
    Largos(i)= DD.wps{i}.length;
    k=k+2;
end

%% Creatting data driven model

Total_volume_ref = sum(model_ref.operators.pv);

L = 3000;
G = cartGrid([10, 1, 36], [L, L/5 ,L/5]*meter^3);
G = computeGeometry(G);

fluid = example.model.fluid;

rock = makeRock(G, 200*milli*darcy, 0.12);

gravity off
model = GenericBlackOilModel(G, rock, fluid);
model.gas=false;
model.OutputStateFunctions = {};
 % Internal faces to be off by setting transmisibility to zero  
 
Total_volume_DD = sum(model.operators.pv)

VolumeRation = Total_volume_DD/Total_volume_ref
 
[model,W,indexs] = createDDmodel(model,10,DD.Graph,W_ref,2);


%% Preparing parameters

WellIP = [];
cell = [];
well_index = []; 
levels = 2;        
for  i = 1:14
    for l = 1 : levels
        well_index = [well_index; i,l];
    end 
    WellIP = [WellIP; W(i).WI(1)];
    WellIP = [WellIP; W(i).WI(2)];  
end


     
Tr_boxlimits = [0.01*TT/2 ; 3*TT/2]';
Pv_boxlimits = [0.01*pv/20 ; 4*pv/20]';
well_boxlimits = [ 0.001*WellIP , ...
                   5*WellIP];
               
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

                                         
          InitalState_cell = struct('name','initSw',...
                                   'type','value',...
                                   'boxLims',[0 1],...
                                   'distribution','cell',...
                                   'Indx',{model.G.cells.indexMap});
                               

parameters =  {};   
parameters{1} = transmisibility_conection;
parameters{2}  = porevolume_conection ;
parameters{3}  = well_IP ;
parameters{4} = InitalState_cell;

% Prepare the model for simulation.

model = model.validateModel();
model.FacilityModel.toleranceWellRate = 1e-4;
model.FacilityModel.toleranceWellMS  = 1e-4;

s0=zeros(360,2);


%I1 -P1
 s0(1:8,1) = 1;
 s0(1:8,2) = 0;
 s0(9:9,1) = 0;
 s0(9:9,2) = 1;
 s0(10,1) = 0;
 s0(10,2) = 1;
  s0(11:17,1) = 1;
 s0(11:17,2) = 0;
 s0(18:20,1) = 0.5;
 s0(18:20,2) = 0.5;

%I2 -P1
lvl = 20;
 s0((1:4)+lvl,1) = 1;
 s0((1:4)+lvl,2) = 0;
 s0((5:10)+lvl,1) = 0;
 s0((5:10)+lvl,2) = 1;
 
 s0((11:15)+lvl,1) = 1;
 s0((11:15)+lvl,2) = 0;
 s0((16:20)+lvl,1) = 0;
 s0((16:20)+lvl,2) = 1;


%I3 -P1
lvl = 40;
 s0((1:5)+lvl,1) = 1;
 s0((1:5)+lvl,2) = 0;
 s0((6:10)+lvl,1) = 0;
 s0((6:10)+lvl,2) = 1;
 
 s0((11:15)+lvl,1) = 1;
 s0((11:15)+lvl,2) = 0;
 s0((16:20)+lvl,1) = 0.5;
 s0((16:20)+lvl,2) = 0.5;

 
 
 %I4 -P1
lvl = 60;
 s0((1:10)+lvl,1) = 0;
 s0((1:10)+lvl,2) = 1;

 s0((11:20)+lvl,1) = 0;
 s0((11:20)+lvl,2) = 1;


  
 %I3 -P2
lvl = 80;
 s0((1:3)+lvl,1) = 1;
 s0((1:3)+lvl,2) = 0;
 s0((4:10)+lvl,1) = 0;
 s0((4:10)+lvl,2) = 1;
 
 s0((11:17)+lvl,1) = 1;
 s0((11:17)+lvl,2) = 0;
 s0((18:20)+lvl,1) = 0.2;
 s0((18:20)+lvl,2) = 0.8;
 
  %I4 -P2
lvl = 100;
 s0((1:10)+lvl,1) = 0;
 s0((1:10)+lvl,2) = 1;

 s0((11:20)+lvl,1) = 0;
 s0((11:20)+lvl,2) = 1;

  %I5 -P2
lvl = 120;
 s0((1:3)+lvl,1) = 1;
 s0((1:3)+lvl,2) = 0;
 s0((4:10)+lvl,1) = 0;
 s0((4:10)+lvl,2) = 1;
 
 s0((11:17)+lvl,1) = 1;
 s0((11:17)+lvl,2) = 0;
 s0((18:20)+lvl,1) = 0;
 s0((18:20)+lvl,2) = 1;
 
 
   %I4 -P3
lvl = 140;
 s0((1:10)+lvl,1) = 0;
 s0((1:10)+lvl,2) = 1;

 s0((11:20)+lvl,1) = 0;
 s0((11:20)+lvl,2) = 1;
 
    %I6 -P3
lvl = 160;
 s0((1:10)+lvl,1) = 0;
 s0((1:10)+lvl,2) = 1;

 s0((11:20)+lvl,1) = 0;
 s0((11:20)+lvl,2) = 1;
 
    %I4 -P4
lvl = 180;
 s0((1:10)+lvl,1) = 0;
 s0((1:10)+lvl,2) = 1;

 s0((11:20)+lvl,1) = 0;
 s0((11:20)+lvl,2) = 1; 
 
 
 
   %I5 -P4
lvl = 200;
 s0((1:3)+lvl,1) = 1;
 s0((1:3)+lvl,2) = 0;
 s0((4:10)+lvl,1) = 0;
 s0((4:10)+lvl,2) = 1;
 
 s0((11:14)+lvl,1) = 1;
 s0((11:14)+lvl,2) = 0;
 s0((15:20)+lvl,1) = 0;
 s0((15:20)+lvl,2) = 1;
 
     %I6 -P4
lvl = 220;
 s0((1:10)+lvl,1) = 0;
 s0((1:10)+lvl,2) = 1;

 s0((11:20)+lvl,1) = 0;
 s0((11:20)+lvl,2) = 1; 
 
 

 
     %I5 -P6

 lvl = 240;
 s0((1:9)+lvl,1) = 1;
 s0((1:9)+lvl,2) = 0;
 s0((10:10)+lvl,1) = 0;
 s0((10:10)+lvl,2) = 1;
 
 s0((11:20)+lvl,1) = 0;
 s0((11:20)+lvl,2) = 1;
 
 
 
     %I6 -P6

 lvl = 260;
 s0((1:10)+lvl,1) = 0;
 s0((1:10)+lvl,2) = 1;
 
 
 s0((11:20)+lvl,1) = 0;
 s0((11:20)+lvl,2) = 1;
 
 %s0((19:20)+lvl,1) = 1;
 %s0((19:20)+lvl,2) = 0;
 
   
 %I7 -P6
lvl = 280;
 s0((1:8)+lvl,1) = 1;
 s0((1:8)+lvl,2) = 0;
 s0((8:10)+lvl,1) = 0.5;
 s0((8:10)+lvl,2) = 0.5;
 
 s0((11:20)+lvl,1) = 1;
 s0((11:20)+lvl,2) = 0;

     %I8 -P6
lvl = 300;
 s0((1:5)+lvl,1) = 1;
 s0((1:5)+lvl,2) = 0;
 s0((6:10)+lvl,1) = 0;
 s0((6:10)+lvl,2) = 1;
 
 s0((11:14)+lvl,1) = 1;
 s0((11:14)+lvl,2) = 0;
 s0((15:20)+lvl,1) = 0;
 s0((15:20)+lvl,2) = 1;
 
 
       %I6 -P5
lvl = 320;
 s0((1:10)+lvl,1) = 0;
 s0((1:10)+lvl,2) = 1;

 s0((11:20)+lvl,1) = 0;
 s0((11:20)+lvl,2) = 1; 
 
    %I8 -P5
lvl = 340;
 s0((1:3)+lvl,1) = 1;
 s0((1:3)+lvl,2) = 0;
 s0((4:10)+lvl,1) = 0;
 s0((4:10)+lvl,2) = 1;
 
 s0((11:14)+lvl,1) = 1;
 s0((11:14)+lvl,2) = 0;
 s0((15:20)+lvl,1) = 0;
 s0((15:20)+lvl,2) = 1;
 
 
 
state0 = initState(model.G, W , 350*barsa,s0); 

figure, plotCellData(model.G, state0.s(:,2)), colorbar, view(0,0);axis equal tight;  daspect([1,0.1,0.1])

dt = schedule_ref.step.val;

 %% Simulation of the reference model and Data-driven model

 schedule = simpleSchedule(dt(1:70), 'W', W);

         weighting =  {'WaterRateWeight',  1e5, ...
                       'OilRateWeight',    1e5, ...
                       'BHPWeight',        1e-3};
 schedule_0=schedule;
                   
 val{1} = TT/2;
 val{2} = pv/20;
 val{3} = WellIP';
 val{4} = state0.s(:,1);
 
  p0_fd = value2control(val,parameters);
   
 [misfitVal_0,gradient,wellSols_0,states_0] = Simulate_BFGS(p0_fd,parameters,model,schedule,state0, wellSols_ref,weighting,1);

          
 plotWellSols({wellSols_ref,wellSols_0},{schedule_ref.step.val,schedule_0.step.val})
          

legend('reference model','initial DD model')



%% 
obj_scaling     = abs(misfitVal_0);      % objective scaling  

objh = @(p)Simulate_BFGS(p,parameters,model,schedule,state0,  wellSols_ref,weighting,obj_scaling);

[v, p_opt, history] = unitBoxBFGS(p0_fd, objh,'gradTol',             1e-2, ...
                                              'objChangeTol',        5e-3);


 schedule = simpleSchedule(dt, 'W', W);



 [misfitVal_opt,gradient_opt,wellSols_opt] = Simulate_BFGS(p_opt,parameters,model,schedule,state0, wellSols_ref,weighting,obj_scaling);
 [misfitVal_0,gradient_0,wellSols_0] = Simulate_BFGS(p0_fd,parameters,model,schedule,state0, wellSols_ref,weighting,obj_scaling);


plotWellSols({wellSols_ref,wellSols_0,wellSols_opt},{schedule_ref.step.val,schedule.step.val,schedule.step.val})
