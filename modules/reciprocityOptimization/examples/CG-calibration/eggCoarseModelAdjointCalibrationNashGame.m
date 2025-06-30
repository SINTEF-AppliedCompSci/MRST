%% Parameter tuning of a very coarse upscaling of the Egg model   
mrstModule add ad-core ad-blackoil deckformat agglom upscaling coarsegrid...
        mrst-gui ad-props incomp optimization test-suite linearsolvers 

%% Setup reference model
test    = TestCase('egg_wo', 'realization', 1);
problem = test.getPackedSimulationProblem();
problem.SimulatorSetup.model.useCNVConvergence = false;
problem.SimulatorSetup.model.nonlinearTolerance=1.0e-9;
problem.SimulatorSetup.model.toleranceMB =1.0e-9;
problem.SimulatorSetup.model.FacilityModel.toleranceWellBHP = 1000;
problem.SimulatorSetup.model.FacilityModel.toleranceWellRate = 5.0e-8;

% Focus on first 60 steps
schedule = problem.SimulatorSetup.schedule;
schedule.step.val     = schedule.step.val(1:60);
schedule.step.control = schedule.step.control(1:60);
Wref     = schedule.control.W;
dt       = schedule.step.val;
nstep    = numel(schedule.step.val);

% Create training schedule with varying controls
perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
rng(0)
schedule = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .01, 'rateFac', .4, 'perturbStep', perturbStep);
problem.SimulatorSetup.schedule = schedule;

% Simulate
simulatePackedProblem(problem);
[wsRef, statesRef] = getPackedSimulatorOutput(problem);
modelRef = problem.SimulatorSetup.model;

% Create schedule with flipped well controls
schedule_new = flipWellControlsAndAddBC(schedule, statesRef, []);
problem_new  = packSimulationProblem(problem.SimulatorSetup.state0, modelRef, schedule_new, 'egg_wo_new');
problem_new.SimulatorSetup.model = modelRef; % Use same model settings
simulatePackedProblem(problem_new);
[wsRef_new, statesRef_new] = getPackedSimulatorOutput(problem_new);
modelRef_new = problem_new.SimulatorSetup.model;

%% Create coarse models for both problems
% Common grid partitioning
blockIx = partitionUI(modelRef.G, [6, 6, 1]);
blockIx = processPartition(modelRef.G, blockIx);
blockIx = compressPartition(blockIx);

% Coarse model for RATE-Inj-BHP-Prod problem
modelCoarse = upscaleModelTPFA(modelRef, blockIx);
modelCoarse.AutoDiffBackend = AutoDiffBackend();
pts = modelCoarse.fluid.krPts;
scaling = {'SWL', pts.w(1), 'SWCR', pts.w(2), 'SWU', pts.w(3), ...
           'SOWCR', pts.ow(2), 'KRW', pts.w(4), 'KRO', pts.ow(4)};
modelCoarse = imposeRelpermScaling(modelCoarse, scaling{:});
modelCoarse.toleranceCNV = 1e-6;

% Coarse model for BHP-Inj-RATE-Prod problem
modelCoarse_new = upscaleModelTPFA(modelRef_new, blockIx);
modelCoarse_new.AutoDiffBackend = AutoDiffBackend();
pts_new = modelCoarse_new.fluid.krPts;
scaling_new = {'SWL', pts_new.w(1), 'SWCR', pts_new.w(2), 'SWU', pts_new.w(3), ...
           'SOWCR', pts_new.ow(2), 'KRW', pts_new.w(4), 'KRO', pts_new.ow(4)};
modelCoarse_new = imposeRelpermScaling(modelCoarse_new, scaling_new{:});
modelCoarse_new.toleranceCNV = 1e-6;

%% Set up initial states and schedules for coarse models
% RATE-Inj-BHP-Prod problem
stateCoarse0 = upscaleState(modelCoarse, modelRef, test.state0);
scheduleCoarse = upscaleSchedule(modelCoarse, schedule, 'wellUpscaleMethod', 'sum');

% BHP-Inj-RATE-Prod problem
stateCoarse0_new = upscaleState(modelCoarse_new, modelRef_new, test.state0);
scheduleCoarse_new = upscaleScheduleNew(modelCoarse_new, schedule_new, modelRef_new, 'wellUpscaleMethod', 'sum');
%% Simulate initial coarse models
[wsCoarse, statesCoarse] = simulateScheduleAD(stateCoarse0, modelCoarse, scheduleCoarse);
[wsCoarse_new, statesCoarse_new] = simulateScheduleAD(stateCoarse0_new, modelCoarse_new, scheduleCoarse_new);
% plot
plotWellSols({wsRef, wsCoarse, wsCoarse_new}, ...
    {schedule.step.val, scheduleCoarse.step.val, scheduleCoarse_new.step.val},...
    'datasetnames',{'fine scale model','initial upscaled RATE-BHP model','initial upscaled BHP-RATE model'});

%% Define parameter configurations for both problems
pv = modelCoarse.operators.pv;
pv_new = modelCoarse_new.operators.pv;

% Parameters for RATE-BHP problem (focus on transmissibility)
config = {...
     %name           include    scaling              boxlims   relativeLimits  
    'porevolume',       1,     'linear',    [.01*pv, 1.5*pv],              []   
    'conntrans',        1,        'log',                  [],     [1e-2, 1e2]      
    'transmissibility', 1,        'log',                  [],     [1e-2, 1e2]  
    'swl',              0,     'linear',             [0, .3],              []
    'swcr',             0,     'linear',             [0, .4],              []
    'swu',              0,     'linear',             [.7, 1],              []
    'sowcr',            0,     'linear',             [0, .4],              []
    'krw',              0,     'linear',           [.5, 1.5],              []
    'kro',              0,     'linear',           [.5, 1.5],              []};

% Parameters for BHP-RATE problem (focus on relperm and pore volume)
config_new = {...
     %name           include    scaling              boxlims   relativeLimits  
    'porevolume',       1,     'linear', [.01*pv_new, 1.5*pv_new],        []   
    'conntrans',        0,        'log',                  [],     [1e-2, 1e2]      
    'transmissibility', 0,        'log',                  [],     [1e-2, 1e2] 
    'swl',              1,     'linear',             [0, .3],              []
    'swcr',             1,     'linear',             [0, .4],              []
    'swu',              1,     'linear',             [.7, 1],              []
    'sowcr',            1,     'linear',             [0, .4],              []
    'krw',              1,     'linear',           [.5, 1.5],              []
    'kro',              1,     'linear',           [.5, 1.5],              []};

%% Initialize parameters for both problems
parameters = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    parameters = addParameter(parameters, struct('model', modelCoarse, 'schedule', scheduleCoarse, 'state0', stateCoarse0), ...
        'name', config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'relativeLimits', config{k,5});
end

parameters_new = [];
for k = 1:size(config_new,1)
    if config_new{k, 2} == 0, continue, end
    parameters_new = addParameter(parameters_new, struct('model', modelCoarse_new, 'schedule', scheduleCoarse_new, 'state0', stateCoarse0_new), ...
        'name', config_new{k,1}, 'scaling', config_new{k,3}, ...
        'boxLims', config_new{k,4}, 'relativeLimits', config_new{k,5});
end

%% Define mismatch function
weighting  = {'WaterRateWeight', 1/(150/day), ...
              'OilRateWeight',   1/(80/day), ...
              'BHPWeight',      1/(20*barsa)};
          
mismatchFn = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref,...
                   'computePartials', compDer, 'tstep', tstep, weighting{:},...
                   'state', state, 'from_states', false);

%% Nash Game Alternating Optimization Algorithm
% Initialize
%% Nash Game Alternating Optimization with Parameter Exchange
% Initialize
setup_init = struct('model', modelCoarse, 'schedule', scheduleCoarse, 'state0', stateCoarse0);
setup_init_new = struct('model', modelCoarse_new, 'schedule', scheduleCoarse_new, 'state0', stateCoarse0_new);

% Get initial parameter vectors
p_opt = getScaledParameterVector(setup_init, parameters);
p_opt_new = getScaledParameterVector(setup_init_new, parameters_new);

% Create combined parameter structure that knows about both problem's parameters
combined_params = struct();
combined_params.rate_bhp = parameters;
combined_params.bhp_rate = parameters_new;

% Optimization settings
maxIt_alt = 5;   % Number of alternations between problems
maxIt = 6;        % Iterations per alternation
tol = 1e-4;       % Tolerance for convergence
history = struct('obj', [], 'obj_new', [], 'params', [], 'params_new', []);
t = 0.5;

% Define objective functions with summands for LM optimization
mismatchFn2 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    matchObservedOW(model, states, schedule, states_ref, ...
        'computePartials', compDer, 'tstep', tstep, weighting{:},...
        'state', state, 'from_states', false, 'mismatchSum', false);

% Main optimization loop
for i = 1:maxIt_alt
    fprintf('\n=== Alternation %d of %d ===\n', i, maxIt_alt);
    
    % 1. Optimize RATE-BHP problem
    fprintf('Optimizing RATE-BHP problem...\n');
    objh2 = @(p) evaluateMatchSummands(p, mismatchFn2, setup_init, combined_params.rate_bhp, statesRef);
    [v2, p_opt, history2] = unitBoxLM(p_opt, objh2, 'maxIt', maxIt);
    
    % Update RATE-BHP model
    setup_init = updateSetupFromScaledParameters(setup_init, combined_params.rate_bhp, p_opt);
    
    % 2. EXCHANGE: Transfer shared parameters to BHP-RATE
    % Get parameter names safely
    rate_bhp_names = arrayfun(@(x) x{1}.name, combined_params.rate_bhp, 'UniformOutput', false);
    bhp_rate_names = arrayfun(@(x) x{1}.name, combined_params.bhp_rate, 'UniformOutput', false);
    
    % Get unscaled parameters
    nparam_rb = cellfun(@(x) x.nParam, combined_params.rate_bhp);
    p_cell_rb = mat2cell(p_opt, nparam_rb, 1);
    
    nparam_br = cellfun(@(x) x.nParam, combined_params.bhp_rate);
    p_cell_br = mat2cell(p_opt_new, nparam_br, 1);
    
    for k = 1:numel(combined_params.rate_bhp)
        current_param = combined_params.rate_bhp{k};
        param_name = current_param.name;
        
        % Find matching parameter in BHP-RATE
        idx = find(strcmp(bhp_rate_names, param_name));
        
        if ~isempty(idx)
            % Unscale and transfer value
            pu_rb = current_param.unscale(p_cell_rb{k});
            pu_br = combined_params.bhp_rate{idx}.unscale(p_cell_br{idx});
            pu = t.*pu_rb + (1-t).*pu_br;
            setup_init_new = combined_params.bhp_rate{idx}.setParameter(setup_init_new, pu);
        end
    end
    
    % 3. Optimize BHP-RATE problem
    fprintf('Optimizing BHP-RATE problem...\n');
    p_opt_new = getScaledParameterVector(setup_init_new, combined_params.bhp_rate);
    objh2_new = @(p) evaluateMatchSummands(p, mismatchFn2, setup_init_new, combined_params.bhp_rate, statesRef_new);
    [v2_new, p_opt_new, history2_new] = unitBoxLM(p_opt_new, objh2_new, 'maxIt', maxIt);
    
    % Update BHP-RATE model
    setup_init_new = updateSetupFromScaledParameters(setup_init_new, combined_params.bhp_rate, p_opt_new);
    
    % 4. EXCHANGE: Transfer back to RATE-BHP
    % Get unscaled parameters
    % Get unscaled parameters
    nparam_rb = cellfun(@(x) x.nParam, combined_params.rate_bhp);
    p_cell_rb = mat2cell(p_opt, nparam_rb, 1);
    
    nparam_br = cellfun(@(x) x.nParam, combined_params.bhp_rate);
    p_cell_br = mat2cell(p_opt_new, nparam_br, 1);
    
    for k = 1:numel(combined_params.bhp_rate)
        current_param = combined_params.bhp_rate{k};
        param_name = current_param.name;
        
        % Find matching parameter in BHP-RATE
        idx = find(strcmp(rate_bhp_names, param_name));
        
        if ~isempty(idx)
            % Unscale and transfer value
            pu_br = current_param.unscale(p_cell_br{k});
            pu_rb = combined_params.rate_bhp{idx}.unscale(p_cell_rb{idx});
            pu = t.*pu_rb + (1-t).*pu_br;
            setup_init = combined_params.rate_bhp{idx}.setParameter(setup_init, pu);
        end
    end
    
    
    % Store history and check convergence
    history.obj = [history.obj history2.val];
    history.obj_new = [history.obj_new history2_new.val];
    history.params = [history.params p_opt'];
    history.params_new = [history.params_new p_opt_new'];
    
    if i > 1 && norm(history.params(:,end) - history.params(:,end-1)) < tol && ...
       norm(history.params_new(:,end) - history.params_new(:,end-1)) < tol
        fprintf('Converged after %d alternations\n', i);
        break;
    end
end

setup_opt2 = updateSetupFromScaledParameters(setup_init, parameters, p_opt);
setup_opt2_new = updateSetupFromScaledParameters(setup_init_new, parameters_new, p_opt_new);

[wellSols_opt2, states_opt2] = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, setup_opt2.schedule);
[wellSols_opt2_new, states_opt2_new] = simulateScheduleAD(setup_opt2_new.state0, setup_opt2_new.model, setup_opt2_new.schedule);

plotWellSols({wsRef_new, wellSols_opt2, wellSols_opt2_new}, ...
repmat({schedule.step.val}, 1, 3), ...
'datasetnames', {'reference','coarse tuned (Q-N)', 'coarse tuned (L-M)'});

%% Compare reference, initial coarse and optimizes coarse for a different schedule 
rng(100);
W = test.schedule.control.W;
s_new = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .01, 'rateFac', .2, 'perturbStep', ones(numel(dt),1));
ws_new1 = simulateScheduleAD(problem.SimulatorSetup.state0, modelRef, s_new, ...
    'NonLinearSolver', problem.SimulatorSetup.NonLinearSolver);

w_opt = setup_opt2.schedule.control(1).W;
w_old = scheduleCoarse.control(1).W;   
s_opt = simpleSchedule(dt, 'W', w_opt);
s_old = simpleSchedule(dt, 'W', w_old);
for kw = 1:numel(Wref)
    s_opt.control.W(kw).val = s_new.control.W(kw).val;
    s_old.control.W(kw).val = s_new.control.W(kw).val;
end
ws_coarse_new1 = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, s_opt);
ws_coarse_old1 = simulateScheduleAD(stateCoarse0, modelCoarse, s_old);
plotWellSols({ws_new1,  ws_coarse_old1, ws_coarse_new1}, ...
    repmat({schedule.step.val}, 1, 3), ...
    'datasetnames', {'reference','coarse initial','coarse tuned'});
