%% Adjoint-Based Parameter Inversion for Egg Model
clear all; close all; clc;
mrstModule add agglom upscaling coarsegrid ad-core ad-blackoil ad-props optimization deckformat test-suite

%% (1) Reference Egg Model Setup
gravity off;
test = TestCase('egg_wo', 'realization', 1);
test.name =  'egg_wo_nogravity';
problemRef = test.getPackedSimulationProblem();
problemRef.BaseName = 'egg_wo_nogravity';
problemRef.OutputHandlers.states.dataFolder = 'egg_wo_nogravity';
problemRef.OutputHandlers.reports.dataFolder = 'egg_wo_nogravity';
problemRef.OutputHandlers.wellSols.dataFolder = 'egg_wo_nogravity';
problemRef.problem.SimulatorSetup.model.gravity = [0 0 0];

% Configure model tolerances
problemRef.SimulatorSetup.model.useCNVConvergence = false;
problemRef.SimulatorSetup.model.nonlinearTolerance = 1.0e-9;
problemRef.SimulatorSetup.model.toleranceMB = 1.0e-9;
problemRef.SimulatorSetup.model.FacilityModel.toleranceWellBHP = 1000;
problemRef.SimulatorSetup.model.FacilityModel.toleranceWellRate = 5.0e-8;
problemRef.SimulatorSetup.model.OutputStateFunctions{end+1} = 'ComponentTotalFlux';

% Focus on first 60 timesteps
scheduleRef = problemRef.SimulatorSetup.schedule;
scheduleRef.step.val = scheduleRef.step.val(1:60);
scheduleRef.step.control = scheduleRef.step.control(1:60);
Wref = scheduleRef.control.W;
dt = scheduleRef.step.val;
nstep = numel(scheduleRef.step.val);

% Create perturbed schedule for training
perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
rng(0); % For reproducibility
scheduleRef = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', 0.01, 'rateFac', 0.4, 'perturbStep', perturbStep);
problemRef.SimulatorSetup.schedule = scheduleRef;

% Simulate reference problem
simulatePackedProblem(problemRef);
[wsRef, statesRef] = getPackedSimulatorOutput(problemRef);
modelRef = problemRef.SimulatorSetup.model;
state0 = problemRef.SimulatorSetup.state0;
%% (2) Duplication Procedure Setup
% This section establishes the dual control systems required for the
% Kohn-Vogelius optimization approach and verifies their equivalence.

% 2.1 Setup boundary conditions for inversion
data = setupReservoirBoundaryConditions(modelRef, statesRef, scheduleRef, 'top', []);

% 2.2 Create original BHP-Rate control system
scheduleBhpRate = scheduleRef;
modelBhpRate = modelRef;
[wellSolsBhpRate, statesBhpRate] = simulateScheduleAD(state0, modelBhpRate, scheduleBhpRate);

% 2.3 Create converted Rate-BHP control system
convMap = struct('bhp', 'rate', 'rate', 'bhp');
scheduleRateBhp = convertWellControls(scheduleRef, statesRef, modelRef, ...
                    'ConversionMap', convMap, 'bc', data.bc);
[wellSolsRateBhp, statesRateBhp] = simulateScheduleAD(state0, modelBhpRate, scheduleRateBhp);

% 2.4 Validate control system equivalence
tolerance = 1e-3;
isValid = true;

for i = 1:numel(wellSolsBhpRate)
    % Compare BHP where applicable
    for w = 1:numel(wellSolsBhpRate{i})
        if isfield(wellSolsBhpRate{i}(w), 'bhp') && isfield(wellSolsRateBhp{i}(w), 'bhp')
            bhpDiff = abs(wellSolsBhpRate{i}(w).bhp - wellSolsRateBhp{i}(w).bhp);
            if bhpDiff > tolerance
                fprintf('BHP mismatch at step %d, well %d: %.2e\n', i, w, bhpDiff);
                isValid = false;
            end
        end
        
        % Compare rates where applicable
        if isfield(wellSolsBhpRate{i}(w), 'qWs') && isfield(wellSolsRateBhp{i}(w), 'qWs')
            rateDiff = abs(wellSolsBhpRate{i}(w).qWs - wellSolsRateBhp{i}(w).qWs);
            if rateDiff > tolerance
                fprintf('Water rate mismatch at step %d, well %d: %.2e\n', i, w, rateDiff);
                isValid = false;
            end
        end
    end
end

if isValid
    fprintf('\nValidation PASSED: Both control systems produce equivalent results within tolerance.\n');
else
    error('Validation FAILED: Control systems produce different results');
end

% 2.5 Visualize control system comparison
summary_plots2 = plotWellSols({wellSolsRef, wellSolsBhpRate}, ...
                            {scheduleRef.step.val, scheduleBhpRate.step.val}, ...
                            'datasetnames', {'Reference Model', 'BHP-Rate Model'});

summary_plots3 = plotWellSols({wellSolsBhpRate, wellSolsRateBhp}, ...
                            {scheduleBhpRate.step.val, scheduleRateBhp.step.val}, ...
                            'datasetnames', {'BHP-Rate Model', 'Rate-BHP Model'});

%% (3) Initial Model Setup for Optimization
% This section prepares the initial guess model that will be used to start
% the optimization procedure.
% 3.0 We construct a prior field
prior.location = data.wellCells;
prior.name = 'porevolume';
prior.value = modelRef.operators.pv(data.wellCells);


% 3.1 Initialize model with estimated porosity
InitRock = modelRef.rock;
initialPoroValue = 0.4;
InitRock.poro = ones(size(InitRock.poro)) * initialPoroValue;

% 3.2 Create initial model with guessed properties
InitModel = GenericBlackOilModel(modelRef.G, InitRock, modelRef.fluid, ...
    'gas', false, 'water', true, 'oil', true);
InitModel.OutputStateFunctions{end+1} = 'ComponentTotalFlux';
InitModel.OutputStateFunctions{end+1} = 'FaceComponentMobility';
InitModel.OutputStateFunctions{end+1} = 'Mobility';

% 3.3 Preserve exact pore volume in well cells from seismic data
seismicdata.cells = data.wellCells;
seismicdata.value = modelRef.operators.pv(seismicdata.cells);
InitModel.operators.pv(seismicdata.cells) = seismicdata.value;

% 3.4 Configure model solver settings for optimization
InitModel.useCNVConvergence = false;
InitModel.nonlinearTolerance = 1.0e-9;
InitModel.FacilityModel.toleranceWellBHP = 1000;
InitModel.FacilityModel.toleranceWellRate = 1.0e-8;
InitModel.FacilityModel.nonlinearTolerance = 1.0e-8;
InitModel.OutputStateFunctions{end+1} = 'ComponentTotalFlux';

% 3.5 Validate initial model with both control systems
modelBhpRate = InitModel;
modelRateBhp = InitModel;

fprintf('\nRunning initial forward simulations with guessed parameters...\n');
[wellSolsBhpRate, statesBhpRate] = simulateScheduleAD(state0, modelBhpRate, scheduleBhpRate);
[wellSolsRateBhp, statesRateBhp] = simulateScheduleAD(state0, modelRateBhp, scheduleRateBhp);

% 3.6 Visualize initial model performance
figure(3);
summary_plots3 = plotWellSols({wellSolsBhpRate, wellSolsRateBhp}, ...
                            {scheduleBhpRate.step.val, scheduleRateBhp.step.val},...
                            'datasetnames', {'BHP-Rate Model (Initial)', ...
                            'Rate-BHP Model (Initial)'});
title('Initial Model Performance Comparison');

%% Optimization Setup Complete
fprintf('\nInitial model setup complete. Ready to begin optimization procedure.\n');
%% Optimization Setup
% Prepare parameter structures for both control systems
params_modelBhpRate = [];  
params_modelRateBhp = [];

% Define model setups
setup_modelBhpRate = struct('model', InitModel, 'schedule', scheduleBhpRate, 'state0', state0);
setup_modelRateBhp = struct('model', InitModel, 'schedule', scheduleRateBhp, 'state0', state0);

% Set parameter bounds based on reference model
pv_range = [min(modelRef.operators.pv).*0+0.1, 0.9+0.*max(modelRef.operators.pv)];
params_modelBhpRate = addParameter(params_modelBhpRate, setup_modelBhpRate, ...
                                 'name', 'porevolume', 'boxLims', pv_range);
params_modelRateBhp = addParameter(params_modelRateBhp, setup_modelRateBhp, ...
                                 'name', 'porevolume', 'boxLims', pv_range);

%% Optimization Configuration
% Define weighting for objective function
weighting = {'WaterRateWeight', (100*stb/day)^-1, ...
             'OilRateWeight',   (100*stb/day)^-1, ...
             'BHPWeight',       (95*barsa)^-1};


%% Two-Stage Optimization Process
% 1. BFGS Optimization (Optional)
use_bfgs = true;  % Set to false to skip BFGS stage
max_iterations = 30;
bfgs_iterations = 10;

if use_bfgs
    fprintf('\nStarting BFGS optimization...\n');
    % First objective function handle
    obj_func_bfgs = @(m1,s1,sch1,m2,s2,sch2,sc,obj_prior, deriv,t,s1_in,s2_in) ...
               matchObservedTwinModelsOW(m1,s1,sch1,m2,s2,sch2,sc,obj_prior,...
               'computePartials',deriv,'tstep',t,'state1',s1_in,...
               'state2',s2_in,'from_states',false,weighting{:});
    p0 = getScaledParameterVector(setup_modelBhpRate, params_modelBhpRate);
    obj_bfgs = @(p) evaluateMatchTwinModels(p, obj_func_bfgs, setup_modelBhpRate, ...
                                          setup_modelRateBhp, params_modelBhpRate, ...
                                          statesRef, prior);
    
    [~, p_opt, ~] = unitBoxBFGS(p0, obj_bfgs, ...
                               'maxIt', bfgs_iterations, ...
                               'objChangeTol', 1e-7, ...
                               'gradTol', 1e-7);
else
    p_opt = getScaledParameterVector(setup_modelBhpRate, params_modelBhpRate);
end

% 2. Levenberg-Marquardt Optimization
obj_func_lm = @(m1,s1,sch1,m2,s2,sch2,sc,obj_prior, deriv,t,s1_in,s2_in) ...
    matchObservedTwinModelsOW(m1,s1,sch1,m2,s2,sch2,sc, obj_prior,...
    'computePartials',deriv,'tstep',t,'state1',s1_in,...
    'state2',s2_in,'from_states',false,weighting{:}, 'mismatchSum', false);
fprintf('\nStarting Levenberg-Marquardt optimization...\n');
obj_lm = @(p) evaluateMatchTwinModelsSummands(p, obj_func_lm, setup_modelBhpRate, ...
                                            setup_modelRateBhp, params_modelBhpRate, ...
                                            statesRef, prior);

[~, p_opt_final, ~] = unitBoxLM(p_opt, obj_lm, 'maxIt', max_iterations);

%% Final Simulation with Optimized Parameters
setup_final = updateSetupFromScaledParameters(setup_modelBhpRate, params_modelBhpRate, p_opt_final);
setup_final.schedule.step = scheduleRef.step; % Use full time horizon

fprintf('\nRunning final simulation with optimized parameters...\n');
[wellSols_opt, states_opt] = simulateScheduleAD(setup_final.state0, ...
                                              setup_final.model, ...
                                              setup_final.schedule);

%% Results Visualization
plotWellSols({wellSolsRef, wellSols_opt}, ...
            {scheduleRef.step.val, setup_final.schedule.step.val}, ...
            'datasetnames', {'Reference Model', 'Optimized Model'});
title('Optimization Results Comparison');

% Create custom "boxplot" using percentiles
ref_pv = modelRef.operators.pv;
init_pv = InitModel.operators.pv;
opt_pv = setup_final.model.operators.pv;
percentiles = [5 25 50 75 95];
ref_prct = prctile(ref_pv, percentiles);
init_prct = prctile(init_pv, percentiles);
opt_prct = prctile(opt_pv, percentiles);

figure;
hold on;
% Reference
plot([1 1], [ref_prct(1) ref_prct(5)], 'b-', 'LineWidth', 2);
rectangle('Position', [0.8, ref_prct(2), 0.4, ref_prct(4)-ref_prct(2)], ...
          'FaceColor', [0.7 0.7 1], 'EdgeColor', 'b');
plot([0.8 1.2], [ref_prct(3) ref_prct(3)], 'b-', 'LineWidth', 2);

% Initial
plot([2 2], [init_prct(1) init_prct(5)], 'r-', 'LineWidth', 2);
rectangle('Position', [1.8, init_prct(2), 0.4, init_prct(4)-init_prct(2)], ...
          'FaceColor', [1 0.7 0.7], 'EdgeColor', 'r');
plot([1.8 2.2], [init_prct(3) init_prct(3)], 'r-', 'LineWidth', 2);

% Optimized
plot([3 3], [opt_prct(1) opt_prct(5)], 'g-', 'LineWidth', 2);
rectangle('Position', [2.8, opt_prct(2), 0.4, opt_prct(4)-opt_prct(2)], ...
          'FaceColor', [0.7 1 0.7], 'EdgeColor', 'g');
plot([2.8 3.2], [opt_prct(3) opt_prct(3)], 'g-', 'LineWidth', 2);

set(gca, 'XTick', 1:3, 'XTickLabel', {'Reference', 'Initial', 'Optimized'});
ylabel('Pore Volume');
title('Distribution Comparison (5-25-50-75-95 percentiles)');
grid on;

% Panel 3: Spatial distribution
figure;
plot(modelRef.operators.pv, 'LineWidth', 2); hold on;
plot(InitModel.operators.pv, '--', 'LineWidth', 1.5);
plot(setup_final.model.operators.pv, 'LineWidth', 1.5);
legend('Reference', 'Initial Guess', 'Optimized', 'Location','best');
xlabel('Cell Index'); ylabel('Pore Volume');
title('Spatial Distribution');
grid on;

% Calculate basic statistics
ref_pv = modelRef.operators.pv;
init_pv = InitModel.operators.pv;
opt_pv = setup_final.model.operators.pv;

stats = struct();
stats.Reference = [mean(ref_pv), std(ref_pv), min(ref_pv), max(ref_pv)];
stats.Initial = [mean(init_pv), std(init_pv), min(init_pv), max(init_pv)];
stats.Optimized = [mean(opt_pv), std(opt_pv), min(opt_pv), max(opt_pv)];

% Display as formatted table
fprintf('%-12s %-8s %-8s %-8s %-8s\n', 'Model', 'Mean', 'Std', 'Min', 'Max');
fprintf('%-12s %-8.2f %-8.2f %-8.2f %-8.2f\n', 'Reference', stats.Reference);
fprintf('%-12s %-8.2f %-8.2f %-8.2f %-8.2f\n', 'Initial', stats.Initial);
fprintf('%-12s %-8.2f %-8.2f %-8.2f %-8.2f\n', 'Optimized', stats.Optimized);

% Plot comparison
figure;
bar([stats.Reference(1); stats.Initial(1); stats.Optimized(1)]);
set(gca, 'XTickLabel', {'Reference', 'Initial', 'Optimized'});
ylabel('Mean Pore Volume');
title('Mean PV Comparison');
grid on;