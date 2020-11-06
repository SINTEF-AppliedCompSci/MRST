mrstModule add ad-core ad-blackoil deckformat diagnostics mrst-gui ad-props incomp optimization

 
%% Run EGG field simulation
%Settup_Egg_simulation 
 


wellSols_ref =  wellSols;
model_ref    = model;
states_ref   = states;
schedule_ref = schedule;
W_ref        = schedule_ref.control.W;

%% Initializing parameters

    pv = Mean_PV;
    TT = Mean_Trans;   

%% Creatting data driven model

L = 435;
G = cartGrid([10, 1, numedges(Graph)], [L, L/5 ,L/5]*meter^3);
G = computeGeometry(G);


fluid =  model_ref.fluid;
rock = makeRock(G, 1000*milli*darcy, 0.1);

gravity off
model = GenericBlackOilModel(G, rock, fluid);
model.gas=false;
model.OutputStateFunctions = {};

 
%[model,W,indexs] = createDDmodel_1(model,10,DD.Graph,W_ref);

obj = GPSNet(model,10,Graph,W_ref);
model = obj.model;
W     = obj.W;
indexs.faces = obj.Graph.Edges.Face_Indices;
indexs.cells = obj.Graph.Edges.Cell_Indices;

figure
subplot(1,2,1)
obj.plotNetwork(model_ref.G,pv)

subplot(1,2,2)
obj.plotNetwork(model_ref.G,TT)


figure
obj.plotGPSNetGrid([])


%% Preparing parameters and scaling values for each one

WellIP = [];
cell = [];
well_index = []; 
levels = 1;        
for  i = 1:numel(W)
    for l = 1 : size(W(i).WI)
        well_index = [well_index; i,l];
        WellIP = [WellIP; 2*W(i).WI(l)];
    end
end


%     pv = 0*Mean_PV +Mean_PV(1) ;
%     TT = 0*Mean_Trans+ Mean_Trans(1);

Tr_boxlimits = [0.01*TT ;1.4*TT]';
Pv_boxlimits = [0.01*pv/10 ; 2*pv/10]';
well_boxlimits = [ 0.01*WellIP , ...
                   7.0*WellIP];
               

               
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
lvl = 0;
s0=zeros( numedges(obj.Graph) ,2);
for i = 1: numedges(obj.Graph)
    nodes= deal(obj.Graph.Edges.EndNodes(i,:));
    
    if (obj.Graph.Nodes.SubWell(nodes(1))==1)
        s0((1:10)+lvl,1) = 0;
        s0((1:10)+lvl,2) = 1;
    else
        s0((1:9)+lvl,1) = 1;
        s0((1:9)+lvl,2) = 0;
        s0((10:10)+lvl,1) = 0;
        s0((10:10)+lvl,2) = 1;
    end
    lvl=lvl+10;
end

model = model.validateModel();
state0 = initState(G, W , 400*barsa,s0); 
dt = schedule.step.val;


 %% Simulating the initial DD model

 schedule = simpleSchedule(dt(1:48), 'W', W);

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

Graph= obj.Graph;

Plot_Parameters