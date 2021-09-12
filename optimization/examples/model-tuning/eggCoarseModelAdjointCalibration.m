%% Parameter tuning of a very coarse upscaling of the Egg model   
mrstModule add ad-core ad-blackoil deckformat ...
               agglom upscaling coarsegrid...
               mrst-gui ad-props incomp optimization...
               example-suite linearsolvers 

%% Setup reference model
example  = MRSTExample('egg_wo', 'realization', 1);
problem  = example.getPackedSimulationProblem();

% We focus on the first 60 steps
schedule = problem.SimulatorSetup.schedule;
schedule.step.val = schedule.step.val(1:60);
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

[wellSolsRef, statesRef] = getPackedSimulatorOutput(problem);
modelRef    = problem.SimulatorSetup.model;



% Plot
%example.plot(statesRef, 'step_index', numel(statesRef))
%plotWellSols(wellSolsRef)
%% Coarse-scale model
% We make a coarse grid defined by a uniform 6 x 6 x 1 partition 
blockIx = partitionUI(modelRef.G, [6, 6, 1]);
blockIx = processPartition(modelRef.G, blockIx);
blockIx = compressPartition(blockIx);
% Perform a simple upscaling to obtain a coarse model
cModel = upscaleModelTPFA(modelRef, blockIx);
cModel.AutoDiffBackend = AutoDiffBackend();
% We want to include rel-perm scaling as tunabale parameters, so include
% these for the coarse model. These parameters have no effect for the
% initial coarse model (they are set equal to the ones given by the
% rel-perm curves).
pts = cModel.fluid.krPts;
scaling = {'SWL',   pts.w(1), 'SWCR', pts.w(2), 'SWU', pts.w(3), ...
           'SOWCR', pts.ow(2), 'KRW',  pts.w(4), 'KRO', pts.ow(4)};
cModel = imposeRelpermScaling(cModel, scaling{:});
cModel.toleranceCNV = 1e-6;  % tighter tolerance to improve gradient accuracy

%% Plot comparison of simulation oputput from the coarse and fine model
figure('position',[100 100 1000 400])
axes('position',[.02 .05 .48 .9]);
plotGrid(modelRef.G, 'EdgeAlpha',.2); 
view(174,60);
title('Fine-scale grid (18553 cells)')
plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10); axis off tight
camlight headlight

axes('position',[.5 .05 .48 .9]);
%plotCellData(cModel.G, cModel.rock.poro, 'EdgeColor', 'none');
plotGrid(cModel.G, 'EdgeAlpha',.8);
title('Coarse-scale grid (33 cells)')
view(174,60); 
plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10); axis off tight
camlight headlight

%% Simulate initial upscaled coarse model for full time
cState0   = upscaleState(cModel, modelRef, example.state0);
cSchedule = upscaleSchedule(cModel, schedule, 'wellUpscaleMethod', 'sum');
[cWellSols, cStates] = simulateScheduleAD(cState0, cModel, cSchedule);

plotWellSols({wellSolsRef, cWellSols}, ...
    {schedule.step.val, cSchedule.step.val},...
    'datasetnames',{'fine scale model','initial upscaled model'}, ...
    'zoom', true, 'field', 'qOs', 'SelectedWells', 7);

%% Specify training schedule and parameters to be matched
% We use the first half of the given data for training. In this setup, we
% use all pore volumes, transmissibilities, and well connections in the
% coarse grid as calibration parameters.
% trainSteps = 1:round(numel(schedule.step.val)/2);
% timeSteps  = schedule.step.val(trainSteps); 
% trainSched = upscaleSchedule(cModel, simpleSchedule(timeSteps, 'W', Wref));
setup_init = struct('model', cModel, 'schedule', cSchedule, 'state0', cState0);

pv = cModel.operators.pv;
% set up 'matrix' for parameter options for easier editing. The specific
% limits set for the various parameters, influences the tuning/optimization 
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
% The mismatch function is defined as a function handle to a library
% function from the optimization module that computes the mismatch between
% a given simulation and a reference state. For an oil-water system, the
% match is computed based on three quantities (water/oil rate and bhp) and
% these should be given an associated weight. Weights on the order the
% reciprocal of rate magnitudes/pressure variations result in a properly 
% scaled mismatch function 
weighting  = {'WaterRateWeight',  1/(150/day), ...
              'OilRateWeight',    1/(80/day), ...
              'BHPWeight',        1/(20*barsa)};
mismatchFn = @(model, states, schedule, states_ref, tt, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
                   'computePartials', tt, 'tstep', tstep, weighting{:},...
                   'state', state, 'from_states', false);

%% Model calibration
pvec = getScaledParameterVector(setup_init, parameters);
objh = @(p) evaluateMatch(p, mismatchFn, setup_init, parameters, statesRef);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
[v, p_opt, history] = unitBoxBFGS(pvec, objh, 'objChangeTol', 1e-5, ...
    'maxIt', 30, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%% Create a new coarse model setup with the optimal parameters, 
%% and rerun the simulation for full time horizon
setup_opt = updateSetupFromScaledParameters(setup_init, parameters, p_opt); 
% reset time-steps to full schedule
% setup_opt.schedule.step = schedule.step;
[wellSols_opt, states_opt] = simulateScheduleAD(setup_opt.state0, setup_opt.model, setup_opt.schedule);

%% Plot well curves for reference fine, initial coarse and tuned coarse models 
fh = plotWellSols({wellSolsRef,cWellSols,wellSols_opt}, ...
    repmat({schedule.step.val}, 1, 3), ...
    'datasetnames', {'reference','coarse initial','coarse tuned'}, 'zoom', true, ...
    'field', 'qOs', 'SelectedWells', 7);

%% Plot the pore volume updates
% fetch pore volume differences in initial and tuned coarse models
dpv = setup_opt.model.operators.pv - setup_init.model.operators.pv;
figure
plotCellData(cModel.G, dpv, 'EdgeColor','none');
plotFaces(cModel.G, boundaryFaces(cModel.G), 'EdgeColor', [0.4 0.4 0.4], ...
         'EdgeAlpha',.5, 'FaceColor', 'none');
view(174,60);
plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10); axis off tight
colorbar('south');
                        
% %% rerun new schedule
% rng(100);
% schedule_new = perturbedSimpleSchedule(dt, 'W', Wref, ...
%     'pressureFac', .01, 'rateFac', .4, 'perturbStep', perturbStep);
% ws_new = simulateScheduleAD(problem.SimulatorSetup.state0, modelRef, schedule_new, ...
%     'NonLinearSolver', problem.SimulatorSetup.NonLinearSolver);
% 
% schedule_new_coarse = setup_opt.schedule;
% schedule_new_old    = cSchedule;
% % set new values
% for kc = 1:numel(schedule_new.control)
%     for kw = 1:numel(Wref)
%         schedule_new_coarse.control(kc).W(kw).val = schedule_new.control(kc).W(kw).val;
%         schedule_new_old.control(kc).W(kw).val = schedule_new.control(kc).W(kw).val;
%     end
% end
% ws_coarse_new = simulateScheduleAD(setup_opt.state0, setup_opt.model, schedule_new_coarse);
% ws_coarse_old = simulateScheduleAD(cState0, cModel, schedule_new_old);
%  plotWellSols({ws_new, ws_coarse_new, ws_coarse_old}, ...
%     repmat({schedule.step.val}, 1, 3), ...
%     'datasetnames', {'reference','tuned', 'original'});
%%
rng(100);
W = example.schedule.control.W;
s_new = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .01, 'rateFac', .2, 'perturbStep', ones(numel(dt),1));
ws_new1 = simulateScheduleAD(problem.SimulatorSetup.state0, modelRef, s_new, ...
    'NonLinearSolver', problem.SimulatorSetup.NonLinearSolver);

w_opt = setup_opt.schedule.control(1).W;
w_old = cSchedule.control(1).W;   
s_opt = simpleSchedule(dt, 'W', w_opt);
s_old = simpleSchedule(dt, 'W', w_old);
for kw = 1:numel(Wref)
    s_opt.control.W(kw).val = s_new.control.W(kw).val;
    s_old.control.W(kw).val = s_new.control.W(kw).val;
end
ws_coarse_new1 = simulateScheduleAD(setup_opt.state0, setup_opt.model, s_opt);
ws_coarse_old1 = simulateScheduleAD(cState0, cModel, s_old);
plotWellSols({ws_new1,  ws_coarse_old1, ws_coarse_new1}, ...
    repmat({schedule.step.val}, 1, 3), ...
    'datasetnames', {'reference','coarse initial','coarse tuned'});

