clearvars;
close all;

%% Initialize MRST Modules
mrstModule add ad-core ad-blackoil deckformat agglom upscaling coarsegrid...
        mrst-gui ad-props incomp optimization test-suite linearsolvers 

%% Setup Reference Model
try
    test = TestCase('egg_wo', 'realization', 1);
    problem = test.getPackedSimulationProblem();
    
    % Configure model tolerances
    problem.SimulatorSetup.model.useCNVConvergence = false;
    problem.SimulatorSetup.model.nonlinearTolerance = 1.0e-9;
    problem.SimulatorSetup.model.toleranceMB = 1.0e-9;
    problem.SimulatorSetup.model.FacilityModel.toleranceWellBHP = 1000;
    problem.SimulatorSetup.model.FacilityModel.toleranceWellRate = 5.0e-8;
    problem.SimulatorSetup.model.OutputStateFunctions{end+1} = 'ComponentTotalFlux';
    
    % Focus on first 60 timesteps
    schedule = problem.SimulatorSetup.schedule;
    schedule.step.val = schedule.step.val(1:60);
    schedule.step.control = schedule.step.control(1:60);
    Wref = schedule.control.W;
    dt = schedule.step.val;
    nstep = numel(schedule.step.val);
    
    % Create perturbed schedule for training
    perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
    rng(0); % For reproducibility
    schedule = perturbedSimpleSchedule(dt, 'W', Wref, ...
        'pressureFac', 0.01, 'rateFac', 0.4, 'perturbStep', perturbStep);
    problem.SimulatorSetup.schedule = schedule;
    
    % Simulate reference problem
    simulatePackedProblem(problem);
    [wsRef, statesRef] = getPackedSimulatorOutput(problem);
    modelRef = problem.SimulatorSetup.model;
    
    % Create schedule with new well controls
    convMap = struct('bhp', 'rate', 'rate', 'bhp');
    schedule_new = convertWellControls(schedule, statesRef, modelRef, 'ConversionMap', convMap);
    problem_new = packSimulationProblem(problem.SimulatorSetup.state0, ...
        problem.SimulatorSetup.model, schedule_new, 'egg_wo_new');
    problem_new.SimulatorSetup.schedule = schedule_new;
    problem_new.Name = 'egg_wo_new';
    simulatePackedProblem(problem_new);
    [wsRef_new, statesRef_new] = getPackedSimulatorOutput(problem_new);
    modelRef_new = problem_new.SimulatorSetup.model;
    
catch ME
    error('Failed during reference model setup: %s', ME.message);
end

%% Create Coarse-Scale Models
try
    % RATE-Inj-BHP-Prod problem
    blockIx = partitionUI(modelRef.G, [12, 12, 1]);
    blockIx = processPartition(modelRef.G, blockIx);
    blockIx = compressPartition(blockIx);
    modelCoarse = upscaleModelTPFA(modelRef, blockIx);
    modelCoarse.AutoDiffBackend = AutoDiffBackend();
    
    % Add rel-perm scaling parameters
    pts = modelCoarse.fluid.krPts;
    scaling = {'SWL', pts.w(1), 'SWCR', pts.w(2), 'SWU', pts.w(3), ...
              'SOWCR', pts.ow(2), 'KRW', pts.w(4), 'KRO', pts.ow(4)};
    modelCoarse = imposeRelpermScaling(modelCoarse, scaling{:});
    
    % BHP-Inj-RATE-Prod problem (use same coarse model)
    modelCoarse_new = modelCoarse;
    
catch ME
    error('Failed during coarse model creation: %s', ME.message);
end

%% Plot Grid Comparison
try
    figure;
    subplot(1,2,1)
    plotCellData(modelRef.G, log10(modelRef.rock.perm(:,1)), 'EdgeAlpha', 0.2); 
    title('Fine-scale grid (18553 cells)')
    plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10);
    axis off tight; view(174,60); camlight headlight;
    
    subplot(1,2,2)
    plotCellData(modelCoarse.G, log10(modelCoarse.rock.perm(:,1)), 'EdgeAlpha', 0.8);
    plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10);
    axis off tight; view(174,60); camlight headlight;
    title('Coarse-scale grid (33 cells)');
catch ME
    warning('Failed during plotting: %s', ME.message);
end

%% Simulate Initial Coarse Models
try
    % RATE-Inj-BHP-Prod problem
    stateCoarse0 = upscaleState(modelCoarse, modelRef, test.state0);
    scheduleCoarse = upscaleSchedule(modelCoarse, schedule, 'wellUpscaleMethod', 'sum');
    [wsCoarse, statesCoarse] = simulateScheduleAD(stateCoarse0, modelCoarse, scheduleCoarse);
    
    % BHP-Inj-RATE-Prod problem
    stateCoarse0_new = upscaleState(modelCoarse_new, modelRef_new, test.state0);
    scheduleCoarse_new = upscaleSchedule(modelCoarse_new, schedule_new, 'wellUpscaleMethod', 'sum');
    [wsCoarse_new, statesCoarse_new] = simulateScheduleAD(stateCoarse0_new, modelCoarse_new, scheduleCoarse_new);
    
    % Plot well solutions
    plotWellSols({wsRef, wsCoarse, wsCoarse_new}, ...
        {schedule.step.val, scheduleCoarse.step.val, scheduleCoarse_new.step.val},...
        'datasetnames', {'fine scale model', 'initial upscaled RATE-BHP model', 'initial upscaled BHP-RATE model'});
catch ME
    error('Failed during coarse model simulation: %s', ME.message);
end

%% Parameter Tuning Setup
try
    % Common setup
    pv = modelCoarse.operators.pv;
    pv_new = modelCoarse_new.operators.pv;
    
    % Parameter configuration
    config = {...
        %name               include  scaling      boxlims            relativeLimits  
        'porevolume',       1,     'linear',    [0.01*pv, 1.5*pv],  []
        'conntrans',        1,     'log',       [],                 [1e-2, 1e2]
        'transmissibility', 1,     'log',       [],                 [1e-2, 1e2]
        'swl',              1,     'linear',    [0, 0.3],           []
        'swcr',             1,     'linear',    [0, 0.4],           []
        'swu',              1,     'linear',    [0.7, 1],           []
        'sowcr',            1,     'linear',    [0, 0.4],           []
        'krw',              1,     'linear',    [0.5, 1.5],         []
        'kro',              1,     'linear',    [0.5, 1.5],         []};
    
    % Create parameter structures
    setup_init = struct('model', modelCoarse, 'schedule', scheduleCoarse, 'state0', stateCoarse0);
    setup_init_new = struct('model', modelCoarse_new, 'schedule', scheduleCoarse_new, 'state0', stateCoarse0_new);
    
    parameters = [];
    for k = 1:size(config,1)
        if config{k, 2} == 0, continue; end
        parameters = addParameter(parameters, setup_init, ...
            'name', config{k,1}, 'scaling', config{k,3}, ...
            'boxLims', config{k,4}, 'relativeLimits', config{k,5});
    end
    
    % Weighting for mismatch function
    weighting = {'WaterRateWeight', 1/(150/day), ...
                'OilRateWeight',   1/(80/day), ...
                'BHPWeight',       1/(20*barsa)};
    
    % Objective functions
    obj = @(model_1, states_1, schedule_1, model_2, states_2, schedule_2, states_c, obj_prior,doPartials, tstep, state1, state2)...
        matchObservedTwinModelsOW(model_1, states_1, schedule_1, model_2, states_2, schedule_2, states_c, obj_prior,...
        'computePartials', doPartials, 'tstep', tstep, 'state1', state1, 'state2', state2, ...
        'from_states', false, weighting{:});
    
    obj2 = @(model_1, states_1, schedule_1, model_2, states_2, schedule_2, states_c, obj_prior, compDer, tstep, state1, state2)...
        matchObservedTwinModelsOW(model_1, states_1, schedule_1, model_2, states_2, schedule_2, states_c, obj_prior, ...
        'computePartials', compDer, 'tstep', tstep, 'state1', state1, 'state2', state2, ...
        'from_states', false, weighting{:}, 'mismatchSum', false);
    
    % Get initial parameter vector
    pvec = getScaledParameterVector(setup_init, parameters);
    objh = @(p)evaluateMatchTwinModels(p, obj, setup_init, setup_init_new, parameters, statesRef, []);
    objh2 = @(p)evaluateMatchTwinModelsSummands(p, obj2, setup_init, setup_init_new, parameters, statesRef, []);
    
catch ME
    error('Failed during parameter setup: %s', ME.message);
end

%% Run Optimization
try
    % Run Levenberg-Marquardt optimization
    [v2_new, p_opt2_new, history2_new] = unitBoxLM(pvec, objh2, 'maxIt', 20);
    
    % Update models with optimized parameters
    setup_opt2 = updateSetupFromScaledParameters(setup_init, parameters, p_opt2_new);
    setup_opt2_new = updateSetupFromScaledParameters(setup_init_new, parameters, p_opt2_new);
    
    % Clear unnecessary fields
    setup_opt2_new.model.FlowDiscretization = [];
    setup_opt2_new.model.FlowPropertyFunctions = [];
    
    % Simulate optimized models
    [wellSols_opt2, states_opt2] = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, setup_opt2.schedule);
    [wellSols_opt2_new, states_opt2_new] = simulateScheduleAD(setup_opt2_new.state0, setup_opt2_new.model, setup_opt2_new.schedule);
    
    % Plot results
    plotWellSols({wsRef_new, wellSols_opt2, wellSols_opt2_new}, ...
        repmat({schedule.step.val}, 1, 3), ...
        'datasetnames', {'reference', 'coarse tuned (Q-N)', 'coarse tuned (L-M)'});
    
    % Plot pore volume updates
    dpv = setup_opt2.model.operators.pv - setup_init.model.operators.pv;
    figure;
    plotCellData(modelCoarse.G, dpv, 'EdgeColor', 'none');
    plotFaces(modelCoarse.G, boundaryFaces(modelCoarse.G), 'EdgeColor', [0.4 0.4 0.4], ...
        'EdgeAlpha', 0.5, 'FaceColor', 'none');
    view(174,60);
    plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10);
    axis off tight;
    colorbar('south');
    
catch ME
    error('Failed during optimization: %s', ME.message);
end

%% Validate on New Schedule
try
    rng(100); % For reproducibility
    W = test.schedule.control.W;
    s_new = perturbedSimpleSchedule(dt, 'W', Wref, ...
        'pressureFac', 0.01, 'rateFac', 0.2, 'perturbStep', ones(numel(dt),1));
    
    % Simulate reference model with new schedule
    [ws_new1, states_new1] = simulateScheduleAD(problem.SimulatorSetup.state0, modelRef, s_new, ...
        'NonLinearSolver', problem.SimulatorSetup.NonLinearSolver);
    
    % Prepare coarse model schedules
    w_opt = setup_opt2.schedule.control(1).W;
    w_old = scheduleCoarse.control(1).W;
    s_opt = simpleSchedule(dt, 'W', w_opt);
    s_old = simpleSchedule(dt, 'W', w_old);
    
    % Apply new control values
    for kw = 1:numel(Wref)
        s_opt.control.W(kw).val = s_new.control.W(kw).val;
        s_old.control.W(kw).val = s_new.control.W(kw).val;
    end
    
    % Simulate coarse models
    ws_coarse_new1 = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, s_opt);
    ws_coarse_old1 = simulateScheduleAD(stateCoarse0, modelCoarse, s_old);
    
    % Plot comparison
    plotWellSols({ws_new1, ws_coarse_old1, ws_coarse_new1}, ...
        repmat({schedule.step.val}, 1, 3), ...
        'datasetnames', {'reference', 'coarse initial', 'coarse tuned'});
    
catch ME
    error('Failed during validation: %s', ME.message);
end