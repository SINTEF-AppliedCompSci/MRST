%% Parameter tuning of a very coarse upscaling of the Egg model   
mrstModule add ad-core ad-blackoil deckformat agglom upscaling coarsegrid...
        mrst-gui ad-props incomp optimization example-suite linearsolvers 

%% Setup reference model
example  = MRSTExample('egg_wo', 'realization', 1);
problem  = example.getPackedSimulationProblem();

% We focus on the first 60 steps
schedule = problem.SimulatorSetup.schedule;
schedule.step.val     = schedule.step.val(1:60);
schedule.step.control = schedule.step.control(1:60);
Wref     = schedule.control.W;
dt       = schedule.step.val;
nstep    = numel(schedule.step.val);

% For the training run, we create a schedule with varying controls
perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
rng(0)
schedule = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .01, 'rateFac', .4, 'perturbStep', perturbStep);
problem.SimulatorSetup.schedule = schedule;

% Simulate
simulatePackedProblem(problem);

[wsRef, statesRef] = getPackedSimulatorOutput(problem);
modelRef    = problem.SimulatorSetup.model;

%% Coarse-scale model
% We make a coarse grid defined by a uniform 6 x 6 x 1 partition 
blockIx = partitionUI(modelRef.G, [6, 6, 1]);
blockIx = processPartition(modelRef.G, blockIx);
blockIx = compressPartition(blockIx);
% Perform a simple upscaling to obtain a coarse model
modelCoarse = upscaleModelTPFA(modelRef, blockIx);
modelCoarse.AutoDiffBackend = AutoDiffBackend();
% We want to include rel-perm scaling as tunabale parameters, so include
% these for the coarse model. These parameters have no effect for the
% initial coarse model (they are set equal to the ones given by the
% rel-perm curves).
pts = modelCoarse.fluid.krPts;
scaling = {'SWL',   pts.w(1), 'SWCR', pts.w(2), 'SWU', pts.w(3), ...
           'SOWCR', pts.ow(2), 'KRW',  pts.w(4), 'KRO', pts.ow(4)};
modelCoarse = imposeRelpermScaling(modelCoarse, scaling{:});
modelCoarse.toleranceCNV = 1e-6;  % tighter tolerance to improve gradient accuracy

%% Plot reference and coarse model grids with wells
figure, subplot(1,2,1)
plotGrid(modelRef.G, 'EdgeAlpha',.2); 
title('Fine-scale grid (18553 cells)')
plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10);
axis off tight, view(174,60), camlight headlight
subplot(1,2,2)
plotGrid(modelCoarse.G, 'EdgeAlpha',.8);
plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10);
axis off tight, view(174,60), camlight headlight
title('Coarse-scale grid (33 cells)')

%% Simulate initial upscaled coarse model and compare to reference
stateCoarse0   = upscaleState(modelCoarse, modelRef, example.state0);
scheduleCoarse = upscaleSchedule(modelCoarse, schedule, 'wellUpscaleMethod', 'sum');
[wsCoarse, statesCoarse] = simulateScheduleAD(stateCoarse0, modelCoarse, scheduleCoarse);

plotWellSols({wsRef, wsCoarse}, ...
    {schedule.step.val, scheduleCoarse.step.val},...
    'datasetnames',{'fine scale model','initial upscaled model'});

%% Specify parameters for tuning
setup_init = struct('model', modelCoarse, 'schedule', scheduleCoarse, 'state0', stateCoarse0);
pv = modelCoarse.operators.pv;
% set up 'matrix' for parameter options for easier editing. The specific
% limits set for the various parameters influences the tuning/optimization 
% procedure to a large extent
config = {...
     %name           include    scaling              boxlims   relativeLimits  
    'porevolume',       1,     'linear',    [.01*pv, 1.5*pv],              []   
    'conntrans',        1,        'log',                  [],     [1e-2, 1e2]      
    'transmissibility', 1,        'log',                  [],     [1e-2  1e2]  
    'swl',              1,     'linear',             [0, .3],              []
    'swcr',             1,     'linear',             [0, .4],              []
    'swu',              1,     'linear',             [.7, 1],              []
    'sowcr',            1,     'linear',             [0, .4],              []
    'krw',              1,     'linear',           [.5, 1.5],              []
    'kro',              1,     'linear',           [.5, 1.5],              []};
parameters = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    parameters = addParameter(parameters, setup_init, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'relativeLimits',config{k,5});
end



%% Define the mismatch function
% Function weighting influences the match of each quantity. Rate-weighting
% should be on the same order as (inverse of) rates. BHP-weighting on the
% order of pressure drop in the model.
weighting  = {'WaterRateWeight',  1/(150/day), ...
              'OilRateWeight',    1/(80/day), ...
              'BHPWeight',        1/(20*barsa)};
% make handle          
mismatchFn = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
                   'computePartials', compDer, 'tstep', tstep, weighting{:},...
                   'state', state, 'from_states', false);

%% Model calibration
pvec = getScaledParameterVector(setup_init, parameters);
objh = @(p) evaluateMatch(p, mismatchFn, setup_init, parameters, statesRef);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
[v, p_opt, history] = unitBoxBFGS(pvec, objh, 'objChangeTol', 1e-5, ...
    'maxIt', 30, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Create a new coarse model setup with the optimized parameters, and rerun 
%  for the optimized parameters
setup_opt = updateSetupFromScaledParameters(setup_init, parameters, p_opt); 
[wellSols_opt, states_opt] = simulateScheduleAD(setup_opt.state0, setup_opt.model, setup_opt.schedule);
% compare reference, initial coarse and optimized coarse model outputs
plotWellSols({wsRef,wsCoarse,wellSols_opt}, ...
              repmat({schedule.step.val}, 1, 3), ...
              'datasetnames', {'reference','coarse initial','coarse tuned'});

%% Plot the pore volume updates
% fetch pore volume differences in initial and tuned coarse models
dpv = setup_opt.model.operators.pv - setup_init.model.operators.pv;
figure
plotCellData(modelCoarse.G, dpv, 'EdgeColor','none');
plotFaces(modelCoarse.G, boundaryFaces(modelCoarse.G), 'EdgeColor', [0.4 0.4 0.4], ...
         'EdgeAlpha',.5, 'FaceColor', 'none');
view(174,60);
plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10); axis off tight
colorbar('south');
                        
%% Compare reference, initial coarse and optimizes coarse for a different schedule 
rng(100);
W = example.schedule.control.W;
s_new = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .01, 'rateFac', .2, 'perturbStep', ones(numel(dt),1));
ws_new1 = simulateScheduleAD(problem.SimulatorSetup.state0, modelRef, s_new, ...
    'NonLinearSolver', problem.SimulatorSetup.NonLinearSolver);

w_opt = setup_opt.schedule.control(1).W;
w_old = scheduleCoarse.control(1).W;   
s_opt = simpleSchedule(dt, 'W', w_opt);
s_old = simpleSchedule(dt, 'W', w_old);
for kw = 1:numel(Wref)
    s_opt.control.W(kw).val = s_new.control.W(kw).val;
    s_old.control.W(kw).val = s_new.control.W(kw).val;
end
ws_coarse_new1 = simulateScheduleAD(setup_opt.state0, setup_opt.model, s_opt);
ws_coarse_old1 = simulateScheduleAD(stateCoarse0, modelCoarse, s_old);
plotWellSols({ws_new1,  ws_coarse_old1, ws_coarse_new1}, ...
    repmat({schedule.step.val}, 1, 3), ...
    'datasetnames', {'reference','coarse initial','coarse tuned'});

