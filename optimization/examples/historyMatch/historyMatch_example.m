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
%   -Fludi parameters related to the scaling of relative permeabilities
%   functions

mrstModule add agglom upscaling coarsegrid...
               ad-core ad-blackoil ad-props...
               optimization deckformat;

use_trans = false;
user_simulation_time = [1 ,1 ,2 5 10 , 10*ones(1,10)]*day(); 
user_training_time_steps = 1:10;

user_partition = [3,3,2];

%% Setting up the fine-scale model
% We make a small model that consists of two different facies with
% contrasting petrophysical properties. 
% Two injectors and two producers are
% placed in the corner at diferent depths.

[model_fine_scale, W, state0] = simpleModelForHMExample();
 model_fine_scale.gas=false;
 model_fine_scale.OutputStateFunctions = {};
 model_fine_scale = model_fine_scale.validateModel();


figure(1), clf,
set(gcf,'Position', [860 450 840 400],'PaperPositionMode','auto');
subplot(1,2,1);
plotCellData(model_fine_scale.G,model_fine_scale.rock.poro,'EdgeAlpha',.5); view(3);
title('Fine-scale model porosity')
plotWell(model_fine_scale.G,W,'Color','k'); axis off tight
pax = [min(model_fine_scale.rock.poro) max(model_fine_scale.rock.poro)];
[hc,hh] = colorbarHist(model_fine_scale.rock.poro, pax,'West');
pos = get(hh,'Position'); set(hh,'Position', pos - [.01 0 0 0]);
pos = get(hc,'Position'); set(hc,'Position', pos - [.03 0 0 0]);
set(gca,'Position',[.12 .075 .315 .9])
drawnow

%% Run simulation for fine-scale model
schedule = simpleSchedule(user_simulation_time, 'W', W);

problem = packSimulationProblem(state0, model_fine_scale, schedule, 'model_fine_scale');
[ok, status] = simulatePackedProblem(problem);
[wellSols_fine_scale, states_fine_scale] = getPackedSimulatorOutput(problem);


%% Coarse-scale model
% We make a 3x3x1 coarse grid. 
    % Permeabilities, pore volumes and well indices are upscaled
    % For the wells, we use specific well conditions and use least squares
    % for the flux.
    
    
G = model_fine_scale.G;
rock =  model_fine_scale.rock;
fluid = model_fine_scale.fluid;

hT = computeTrans(G, rock);

p  = partitionCartGrid(G.cartDims, user_partition );
CG = coarsenGeometry(generateCoarseGrid(G, p));

crock = convertRock2Coarse(G, CG, rock);
crock.perm = upscalePerm(G, CG, rock);
[~,~,WC]   = upscaleTrans(CG, hT, 'match_method', 'lsq_flux', ...
                          'bc_method', 'wells_simple', 'wells', W);    
                      


figure(1)
subplot(1,2,2); cla;

plotCellData(CG,crock.poro,'EdgeColor','none');
title('Coarse-scale model porosity')

plotFaces(CG, boundaryFaces(CG),...
   'EdgeColor', [0.4 0.4 0.4],'EdgeAlpha',.5, 'FaceColor', 'none'); view(3);
plotWell(G,W,'Color','k'); axis off tight
[hc,hh]=colorbarHist(crock.poro, pax,'East');
pos = get(hh,'Position'); set(hh,'Position', pos + [.01 0 0 0]);
pos = get(hc,'Position'); set(hc,'Position', pos + [.02 0 0 0]);
set(gca,'Position',[.56 .075 .315 .9])
drawnow

 %% Simulating upscaled coarse-scale model
 % After coarsening we create a coarse model and we compare the production
 % results between the fine scaled and coarse scale models
 
fluid_2 = fluid;
fluid_2 .krPts  = struct('w', [0 0 1 1], 'ow', [0 0 1 1]);
scaling = {'SWL', .4, 'SWCR', .4, 'SWU', .8, 'SOWCR', .2, 'KRW', .6, 'KRO', .6};


model_coarse_scale = GenericBlackOilModel(CG, crock, fluid_2);
model_coarse_scale.gas=false;

if use_trans
    %model_coarse_scale= upscaleModelTPFA(model_fine_scale, user_partition, 'transCoarse',CTrans(model_coarse_scale.operators.internalConn));
    [~,CTrans] = upscaleTrans(CG, hT, 'match_method', 'max_flux', ...
                          'bc_method', 'bc_simple');
    model_coarse_scale.operators.T = CTrans(model_coarse_scale.operators.internalConn);
    model_coarse_scale.AutoDiffBackend = AutoDiffBackend;
end

% Expanding for fluid parameters calibration. 
model_coarse_scale = imposeRelpermScaling(model_coarse_scale, scaling{:});
model_coarse_scale.OutputStateFunctions = {};
model_coarse_scale = model_coarse_scale.validateModel();
model_coarse_scale.toleranceCNV = 1e-6;

state0_coarse = upscaleState(model_coarse_scale, model_fine_scale, state0);

% schedule_training = upscaleSchedule(model_coarse_scale, schedule);
schedule_training = simpleSchedule(user_simulation_time(user_training_time_steps), 'W', WC);

[wellSols_coarse_scale,states_coarse_scale] = simulateScheduleAD(state0_coarse, model_coarse_scale, schedule_training);

summary_plots = plotWellSols({wellSols_fine_scale,wellSols_coarse_scale},{schedule.step.val,schedule_training.step.val});
drawnow

%% Preparing parameters and scaling values for each one
% We define each parameters it's scailing values.

SimulatorSetup = struct('model', model_coarse_scale, 'schedule', schedule_training, 'state0', state0_coarse);

n_cells =  model_coarse_scale.G.cells.num;
% Fluid Parameters
 parameters{1} = ModelParameter(SimulatorSetup, 'name', 'swl','lumping',ones(n_cells,1),'boxLims',[0.00 0.5]);
 parameters{2} = ModelParameter(SimulatorSetup, 'name', 'swcr','lumping',ones(n_cells,1),'boxLims',[0.0 0.5]);
 parameters{3} = ModelParameter(SimulatorSetup, 'name', 'swu','lumping',ones(n_cells,1),'boxLims',[0.55 1.0]);
 parameters{4} = ModelParameter(SimulatorSetup, 'name', 'kro','lumping',ones(n_cells,1),'boxLims',[0.6 1.0]);
 parameters{5} = ModelParameter(SimulatorSetup, 'name', 'krw','lumping',ones(n_cells,1),'boxLims',[0.6 1.0]);

 % Well, porevolume and transmisibility
 parameters{6} = ModelParameter(SimulatorSetup, 'name', 'conntrans','relativeLimits', [.01 4]);
 parameters{7} = ModelParameter(SimulatorSetup, 'name', 'porevolume','relativeLimits', [.01 4]);
 parameters{8} = ModelParameter(SimulatorSetup, 'name', 'transmissibility','relativeLimits', [.01 4]);


 %% Optimization :  History Matching

values = applyFunction(@(p)p.getParameterValue(SimulatorSetup), parameters);
% scale values
u = cell(size(values));
for k = 1:numel(u)
    u{k} = parameters{k}.scale(values{k});
end
p0_ups = vertcat(u{:});  
  
% Defining the weights to evaluate the match
weighting =  {'WaterRateWeight',  (5/day)^-1, ...
              'OilRateWeight',    (5/day)^-1,...
              'BHPWeight',        (50*barsa)^-1};            

 obj = @(model, states, schedule_training, states_ref, tt, tstep, state) matchObservedOW(model, states, schedule_training, states_fine_scale,...
           'computePartials', tt, 'tstep', tstep, weighting{:},'state',state,'from_states',false);
 
objScaling     =  1; % objective scaling  
objh = @(p)evaluateMatch(p, obj, state0_coarse, model_coarse_scale, schedule_training, objScaling ,parameters,  states_fine_scale );

clf(figure(10),'reset');
[v, p_opt, history] = unitBoxBFGS(p0_ups, objh,'objChangeTol',  1e-8, 'maxIt', 25, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Running the simulations to compare initial model vs calibrated model
schedule = simpleSchedule(user_simulation_time, 'W', WC);
[misfitVal_opt,~,wellSols_opt]    = evaluateMatch(p_opt, obj, state0_coarse, model_coarse_scale, schedule, objScaling ,parameters, states_fine_scale,'Gradient','none');
[misfitVal_0  ,~,wellSols_0]      = evaluateMatch(p0_ups,obj, state0_coarse, model_coarse_scale, schedule, objScaling, parameters, states_fine_scale,'Gradient','none');                                                          

try
    figure(summary_plots.Number)
catch
    summary_plots = figure;
end
plotWellSols({wellSols_fine_scale,wellSols_0,wellSols_opt},...
              {schedule.step.val,schedule.step.val,schedule.step.val},...
              'datasetnames',{'fine scale model','initial upscaled model','history matched upscaled model'},...
              'linestyles',{'o', '--', '-'},...
              'figure',summary_plots.Number)
drawnow

