mrstModule add ad-core ad-blackoil deckformat diagnostics mrst-gui ad-props incomp optimization linearsolvers

 
%
load 'Brugge_Model_run_sintef.mat'


wellSols_ref =  wellSols;
model_ref    = model;
states_ref   = {states{75:end}}';
state0_ref   = initState(model.G, schedule.control.W , states{74}.pressure,states{74}.s); 
schedule_ref = schedule;
schedule_ref.step.val= schedule_ref.step.val(75:end);
schedule_ref.step.control= schedule_ref.step.control(75:end);

W_ref        = schedule.control.W;


%  
DD = WellPairNetwork(model_ref,schedule_ref,states_ref,state0_ref,wellSols_ref);

 DD =  DD.filter_wps(150*stb/day);
       DD.plotWellPairConnections()
       DD.plotWellPairsData('subplot',[7,10])       

 % Initializing parameters transmisibility and pore volume, 2 for each well
 % pair
       
k =1;
for i =  1:numel(DD.wps)
    pv([k k+1]) = DD.wps{i}.volume;
    TT([k k+1]) = DD.wps{i}.Tr;
    Largos(i)= DD.wps{i}.length;
    k=k+2;
end



%% Creatting data driven model

L = 2750;
G = cartGrid([10, 1, 2*numel(DD.wps)], [L, L/5 ,L/5]*meter^3);
G = computeGeometry(G);

%fluid = model_ref.fluid;
fluid = initSimpleADIFluid('mu',    [1, 5]*centi*poise, ...    % Transformd into Corey function
                           'rho',   [1000, 900]*kilogram/meter^3, ...
                           'n',     [2, 2], ...
                           'cR',    1e-8/barsa, ...
                           'phases', 'wo');
%CoreyFluid is giving  convergence problem in CNV                     
% fluid_2 = initCoreyFluid('mu',    [0.3200, 1.2940]*centi*poise, ...    % Transformd into Corey function
%                            'rho',   [ 1.0033e+03, 896.9539]*kilogram/meter^3, ...
%                            'n',     [2, 3], ...
%                            'sr',  [0.2221 0 ],...
%                            'kwm',  [1, 0.4]);
                                              
                       
fluid.bO = @(p, varargin) exp((p/barsa - 100)*0.001);

rock = makeRock(G, 200*milli*darcy, 0.2);

gravity off
%model = GenericBlackOilModel(G, rock, fluid);
model = TwoPhaseOilWaterModel(G, rock, fluid);
%model.gas=false;
model.OutputStateFunctions = {};

[model,W,indexs] = createDDmodel(model,10,DD.Graph,W_ref,2);


%% Preparing parameters and scaling values for each one

WellIP = [];
cell = [];
well_index = []; 
levels = 2;        
for  i = 1:numel(W)
    for l = 1 : levels
        well_index = [well_index; i,l];
    end
 
    WellIP = [WellIP; W(i).WI(1)];
    WellIP = [WellIP; W(i).WI(2)];  
end

% scaling parameters
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
%parameters{4} = InitalState_cell;


% Prepare the model for simulation.

model = model.validateModel();
%model.FacilityModel.toleranceWellRate = 1e-3;
%model.FacilityModel.toleranceWellMS  = 1e-3;
model.toleranceCNV = 1e-8;

% initializing saturation manually
InitialSaturation_75

%% Prepare the model for simulation.


state0 = initState(model.G, W , 155*barsa,s0); 
figure, plotCellData(model.G, state0.s(:,2)), colorbar, view(0,0);axis equal tight;  daspect([1,0.1,0.1])

dt = schedule_ref.step.val;

 %% Simulating the initial DD model  in the training period

 schedule = simpleSchedule(dt(1:20), 'W', W);


         weighting =  {'WaterRateWeight',  1e2, ...
                       'OilRateWeight',    2*1e2, ...
                       'BHPWeight',        1e-6};
 schedule_0=schedule;
 
      
 val{1} = TT/2;
 val{2} = pv/20;
 val{3} = WellIP';
 %val{4} = state0.s(:,1);
 
  p0_fd = value2control(val,parameters);
  
  [model,schedule,state0] = control2problem(p0_fd,model,schedule,state0, parameters);

  
 [misfitVal_0,gradient,wellSols_0,states_0] = Simulate_BFGS(p0_fd,parameters,model,schedule,state0, wellSols_ref,weighting,1);

          
 plotWellSols({wellSols_ref,wellSols_0},{schedule_ref.step.val,schedule_0.step.val})
 

  %% Optimization in the training period

obj_scaling     = abs(misfitVal_0);      % objective scaling  



objh = @(p)Simulate_BFGS(p,parameters,model,schedule,state0,  wellSols_ref,weighting,obj_scaling);

[v, p_opt, history] = unitBoxBFGS(p0_fd, objh,'gradTol',             1e-2, ...
                                              'objChangeTol',        5e-3);


 schedule = simpleSchedule(dt, 'W', W);


 [misfitVal_opt,gradient_opt,wellSols_opt] = Simulate_BFGS(p_opt,parameters,model,schedule,state0, wellSols_ref,weighting,obj_scaling);
 [misfitVal_0,gradient_0,wellSols_0] = Simulate_BFGS(p0_fd,parameters,model,schedule,state0, wellSols_ref,weighting,obj_scaling);
 

plotWellSols({wellSols_ref,wellSols_0,wellSols_opt},{schedule_ref.step.val,schedule.step.val,schedule.step.val})
legend('reference model','initial DD model','optimize DD model')
%legend('reference model','initial DD model','optimize DD model')

 %% 
 PlotBrugge_results