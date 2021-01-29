%% Illustrate use Adjoin optimization for history match
% In this example we explain how to history match two cosared models. One 
% model that have been obtained by upscailing an fine scaled model. This
% upscaled coarse model has production rates very similar to the one
% obtained by the fine scale model however the diference can be minimize by
% a history match. 

% A second coarse model with production result diferent from the fine
% scaled model is also used to show how the robust is the adjoint method
% for history matching

% We history match water and oil production rates and we define as parameters:
%   -Transmisibility of each internal faces, 
%   -Pore volume of each cell
%   -Well index for each well 

mrstModule add agglom upscaling coarsegrid...
               ad-core ad-blackoil ad-props...
               optimization;

use_trans = true;

%% Make model
% We make a small model that consists of two different facies with
% contrasting petrophysical properties. 
% Two injectors and two producers are
% placed in the corner at diferent depths.
G  = computeGeometry(cartGrid([40 20 15]));
K1 = gaussianField(G.cartDims, [200 2000]);
p1 = K1(:)*1e-4 + .2;
K1 = K1(:)*milli*darcy;
K2 = gaussianField(G.cartDims, [10 500]);
p2 = K2(:)*1e-4 + .2;
K2 = K2(:)*milli*darcy;

rad1 = G.cells.centroids(:,1).^2 + .5*G.cells.centroids(:,2).^2 ...
   + (G.cells.centroids(:,3)-2).^2;
rad2 = .5*(G.cells.centroids(:,1)-40).^2 + 2*G.cells.centroids(:,2).^2 ...
   + 2*(G.cells.centroids(:,3)-2).^2;

ind = ((rad1>600) & (rad1<1500)) | ((rad2>700) & (rad2<1400));
rock.perm = K2(:);
rock.perm(ind) = K1(ind);
rock.perm = bsxfun(@times, rock.perm, [1 1 1]);
rock.poro = p2(:);
rock.poro(ind) = p1(ind);
pv = poreVolume(G, rock);

W = verticalWell([], G, rock, 4, 17, 4:10, 'Comp_i', [1 0], ...
   'Type', 'rate', 'Val', 100*stb/day, 'Name', 'I1');
W = verticalWell(W, G, rock, 2, 3, 1:6, 'Comp_i', [1 0], ...
    'Type', 'rate', 'Val', 100*stb/day, 'Name', 'I2');
W = verticalWell(W,  G, rock, 35, 3, 1:8, 'Comp_i', [0 1], ...
   'Type', 'bhp', 'Val', 95*barsa, 'Name', 'P1');
W = verticalWell(W,  G, rock, 35, 17, 8:15, 'Comp_i', [0 1], ...
   'Type', 'bhp', 'Val', 95*barsa, 'Name', 'P2');

figure(1); clf,
set(gcf,'Position', [860 450 840 400],'PaperPositionMode','auto');
subplot(1,2,1);
plotCellData(G,rock.poro,'EdgeAlpha',.5); view(3);
title('Porosity')
plotWell(G,W,'Color','k'); axis off tight
pax = [min(rock.poro) max(rock.poro)];
[hc,hh] = colorbarHist(rock.poro, pax,'West');
pos = get(hh,'Position'); set(hh,'Position', pos - [.01 0 0 0]);
pos = get(hc,'Position'); set(hc,'Position', pos - [.03 0 0 0]);
set(gca,'Position',[.12 .075 .315 .9])

%% Setting up the fine-scale model
% Transmissibility
hT = computeTrans(G, rock);
trans = 1 ./ accumarray(G.cells.faces(:,1), 1 ./ hT, [G.faces.num, 1]);

% Fluid model
fluid = initSimpleADIFluid('phases', 'WO',... %Fluid phase
                           'mu' , [1, 10]*centi*poise     , ...%Viscosity
                           'rho', [1014, 859]*kilogram/meter^3, ...%Surface density [kg/m^3]
                           'n', [2 2]);
% Optional fluid parameters                                              
                           %'c', 1e-4/barsa,... %Fluid compressibility
                           %'cR', 1e-5/barsa); % Rock compressibility
                           
                           %fluid.bO = @(p, varargin) exp((p/barsa - 100)*0.01);
                           
% Initial state
s0 = [0.2,0.8];
p0 = 100*barsa;
state0  = initState(G, W,p0,s0);

% Create fine scaled model
model_fine_scale = GenericBlackOilModel(G, rock, fluid);
model_fine_scale.gas=false;
model_fine_scale.OutputStateFunctions = {};
model_fine_scale = model_fine_scale.validateModel();

 %% Simulating on fine-scale model
 
dt = [1 ,1 ,2 5 10 , 10*ones(1,10)]*day();
schedule = simpleSchedule(dt, 'W', W);


%% Run simulation
problem = packSimulationProblem(state0, model_fine_scale, schedule, 'model_fine_scale');
[ok, status] = simulatePackedProblem(problem);
[wellSols_fine_scale, states_fine_scale] = getPackedSimulatorOutput(problem);


%% Coarse-scale model
% We make a 3x3x1 coarse grid. 
% TODO: We are upscailing permeability make the same with transmisibility 
% and update this comments:

    % Transmissibilities and well indices are upscaled
    % using two slightly different methods: for the transmissibilities we use
    % global generic boundary conditions and on each coarse face use the
    % solution that has the largest flux orthogonal to the face to compute the
    % upscaled transmissibility. For the wells, we use specific well
    % conditions and use least squares for the flux.
p  = partitionCartGrid(G.cartDims, [3 3 1]);
CG = coarsenGeometry(generateCoarseGrid(G, p));

crock = convertRock2Coarse(G, CG, rock);
crock.perm = upscalePerm(G, CG, rock);
[~,~,WC]   = upscaleTrans(CG, hT, 'match_method', 'lsq_flux', ...
                          'bc_method', 'wells_simple', 'wells', W);                     

figure(1); subplot(1,2,2); cla
plotCellData(CG,crock.poro,'EdgeColor','none');
plotFaces(CG, boundaryFaces(CG),...
   'EdgeColor', [0.4 0.4 0.4],'EdgeAlpha',.5, 'FaceColor', 'none'); view(3);
plotWell(G,W,'Color','k'); axis off tight
[hc,hh]=colorbarHist(crock.poro, pax,'East');
pos = get(hh,'Position'); set(hh,'Position', pos + [.01 0 0 0]);
pos = get(hc,'Position'); set(hc,'Position', pos + [.02 0 0 0]);
set(gca,'Position',[.56 .075 .315 .9])


 %% Simulating on coarse-scale model
 % After upscailing we create a coarse model and we compare the production
 % results between the fine scaled and coarse scale models
model_coarse_scale = GenericBlackOilModel(CG, crock, fluid);
model_coarse_scale.gas=false;
model_coarse_scale.OutputStateFunctions = {};
model_coarse_scale = model_coarse_scale.validateModel();
 
state0  = initState(CG, WC,p0,s0);

schedule = simpleSchedule(dt, 'W', WC);
[wellSols_coarse_scale,states_coarse_scale] = simulateScheduleAD(state0, model_coarse_scale, schedule);

figure
plotWellSols({wellSols_fine_scale,wellSols_coarse_scale},{schedule.step.val,schedule.step.val})

%% Preparing parameters and scaling values for each one
% We define each parameters well index transmisibility, pore volume, and
% it's scailing values.


% Well index
WellIP = [];
cell = [];
well_index = []; 
levels = 1;
well_index_2 = {};
for  i = 1:numel(WC)
    for l = 1 : levels
        well_index = [well_index; i,l];
        well_index_2{i} = i;
    end
    
    WellIP = [WellIP; WC(i).WI(1)];
end

% Transmisibility and porosity 
TT =  model_coarse_scale.operators.T;
pv =  model_coarse_scale.operators.pv;

Tr_boxlimits = [0.1*min(TT) ; 1.2*max(TT)]';
Pv_boxlimits = [0.1*min(pv); 1.2*max(pv)]';
well_boxlimits = [ 0.1*WellIP , ...
                     2*WellIP];

              
          well_IP = struct('name','conntrans',...
                                   'type','value',...
                                   'boxLims', well_boxlimits,...
                                   'distribution','general',...
                                   'Indx',well_index);


          transmisibility_conection = struct('name','transmissibility',...
                                   'type','value',...
                                   'boxLims', Tr_boxlimits ,...
                                   'distribution','cell',...
                                   'Indx',1:numel(TT));
                               

          porevolume_conection = struct('name','porevolume',...
                                   'type','value',...
                                   'boxLims',Pv_boxlimits,...
                                   'distribution','cell',...
                                   'Indx',1:numel(pv));                                                                
                               
parameters =  {};                          
parameters{1} = transmisibility_conection;
parameters{2}  = porevolume_conection ;
parameters{3}  = well_IP ;

 %% Optimization 1: upscaled coarse model
 
 
% Scailing the parameters to the interval [0,1]
  val = {};              
  val{1} = TT;
  val{2} = pv;
  val{3} = WellIP;
  p0_ups = value2control(val,parameters);
 
% Defining the weights to evaluate the match
  weighting =  {'WaterRateWeight',  1e6, ...
                'OilRateWeight',    1e6, ...
                'BHPWeight',        1e-5};

 obj = @(model, states, schedule, states_ref, tt, tstep, state) matchObservedOW(model, states, schedule, states_fine_scale,...
           'computePartials', tt, 'tstep', tstep, weighting{:},'state',state,'from_states',false);

 objScaling = 1;       
 [misfitVal_0,gradient,wellSols_0,states_0] = evaluateMatch(p0_ups,obj,state0,model_coarse_scale,schedule,objScaling,parameters, states_fine_scale);

                                                           

  
obj_scaling     = abs(misfitVal_0);      % objective scaling  
objh = @(p)evaluateMatch(p, obj, state0, model_coarse_scale, schedule, obj_scaling ,parameters,  states_fine_scale);

[v, p_opt, history] = unitBoxBFGS(p0_ups, objh,'gradTol',             1e-2, ...
                                              'objChangeTol',        0.5e-3,...
                                              'wolfe1',             0.5e-0, ...
                                              'wolfe2',              0.9);
[misfitVal_opt,gradient_opt,wellSols_opt] = evaluateMatch(p_opt, obj, state0, model_coarse_scale, schedule, obj_scaling ,parameters, states_fine_scale);
 

plotWellSols({wellSols_fine_scale,wellSols_0,wellSols_opt},...
              {schedule.step.val,schedule.step.val,schedule.step.val},...
              'datasetnames',{'fine scale model','initial upscaled model','history matched upscaled model'},...
              'linestyles', {'o', '--', '-'})

 %% Optimization 2:  initial upscaled model with a ramdom perturbation
 
 
% Scailing the parameters to the interval [0,1]
  val = {};              
  val{1} = TT+(rand(size(TT))-0.5).*TT;
  val{2} = pv+(rand(size(pv))-0.5).*pv;
  val{3} = WellIP+ (rand(size(WellIP))-0.5).*WellIP;
  p0_mean = value2control(val,parameters);
 
       
 [misfitVal_0,gradient,wellSols_0_mean,states_0] = evaluateMatch(p0_mean,obj,state0,model_coarse_scale,schedule,objScaling,parameters, states_fine_scale);

obj_scaling     = abs(misfitVal_0);      % objective scaling  
objh = @(p)evaluateMatch(p,obj, state0, model_coarse_scale, schedule, obj_scaling ,parameters,  states_fine_scale);

[v, p_opt, history] = unitBoxBFGS(p0_mean, objh,'gradTol',             1e-2, ...
                                              'objChangeTol',        0.5e-3,...
                                              'wolfe1',             0.5e-0, ...
                                              'wolfe2',              0.9);
[misfitVal_opt,gradient_opt,wellSols_opt_mean] = evaluateMatch(p_opt, obj, state0, model_coarse_scale, schedule, obj_scaling ,parameters, states_fine_scale);


plotWellSols({wellSols_fine_scale,wellSols_0,wellSols_opt,wellSols_0_mean,wellSols_opt_mean},...
              {schedule.step.val,schedule.step.val,schedule.step.val,schedule.step.val,schedule.step.val},...
              'datasetnames',{'fine scale model','initial upscaled model','history matched upscaled model','perturbed upscaled model','history matched perturbed upscaled model'},...
              'linestyles', {'o', '--', '-','--', '-'})
