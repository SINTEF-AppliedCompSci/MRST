%% Comprehensive Hydrogen Storage: Black-Oil vs Compositional Comparison
% ===========================================================================
% This script compares three modeling approaches using the same dome-shaped
% aquifer setup from the black-oil example:
% 1. Black-oil model with tabulated PVT data (original script)
% 2. Compositional model without bacterial effects
% 3. Compositional model with bacterial growth and bio-clogging
%
% The simulation uses identical grid, rock properties, wells, and schedule
% for fair comparison between different physical models.
% ---------------------------------------------------------------------------

clearvars; clc;
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui deckformat

%% ============ PART 1: BLACK-OIL MODEL (Original Script) =====================
fprintf('=== Setting up Black-Oil Model ===\n');
simpleBlackOil2DModel();

%% ============ PART 2: COMPOSITIONAL MODELS =====================
fprintf('=== Setting up Compositional Models ===\n');

%% Define Compositional Fluid (H2, H2O, CO2, CH4)
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, ...
    {'H2O', 'H2', 'CO2', 'C1'});

%% Common Parameters for Compositional Models
z0 = [0.995, 0.0, 0.005, 0.0];  % Initial composition [H2O, H2, CO2, C1]
T0 =  273.15 +60;
%% Model 1: Compositional without Bacteria
backend = DiagonalAutoDiffBackend('modifyOperators', true);
model_comp_nobact = BiochemistryModel(G, rock, fluid_bo, compFluid, ...
    true, backend, 'water', false, 'oil', true, 'gas', true, ...
    'bacteriamodel', false, 'liquidPhase', 'O', 'vaporPhase', 'G');
model_comp_nobact.gravity = grav;

% Initialize compositional state without bacteria
state0_comp_nobact = initCompositionalState(model_comp_nobact, pInit, T0, [1.0, 0.0], z0);

%% Model 2: Compositional with Bacteria (if enabled)
model_comp_bact = BiochemistryModel(G, rock, fluid_bo, compFluid, ...
    true, backend, 'water', false, 'oil', true, 'gas', true, ...
    'bacteriamodel', true, 'liquidPhase', 'O', 'vaporPhase', 'G');
model_comp_bact.gravity = grav;
nbact0 = 1;
model_comp_bact.nbactMax = 1e8;
% Initialize with bacteria
state0_comp_bact = initCompositionalStateBacteria(model_comp_bact, pInit, T0, [1.0, 0.0], z0, nbact0);

%% Adjust Schedule for Compositional Models (Well Components)
schedule_comp = schedule_bo;
schedule_comp.control(1).W.components = [0.0, 0.95, 0.05, 0.0];
schedule_comp.control(1).W.T = T0;
schedule_comp.control(1).bc.components = [1.0, 0.0, 0.0, 0.0];

%% ============ PART 4: RUN SIMULATIONS =====================
fprintf('=== Running Simulations ===\n');

%% Setup Nonlinear Solver (Common for all models)
nls = NonLinearSolver();
lsolve = selectLinearSolverAD(model_bo);
nls.LinearSolver = lsolve;

%% Run Black-Oil Simulation
fprintf('1. Running Black-Oil simulation...\n');
clf;
fn = getPlotAfterStep(state0_bo, model_bo, schedule_bo, ...
    'view', [0, 180], 'field', 's:2');
[ws_bo, states_bo, report_bo] = simulateScheduleAD(state0_bo, model_bo, schedule_bo, 'afterStepFn', fn);

%% Run Compositional without Bacteria
fprintf('2. Running Compositional simulation (no bacteria)...\n');
problem_comp_nobact = packSimulationProblem(state0_comp_nobact, model_comp_nobact, schedule_comp, ...
    'H2Storage_Comp_NoBact', 'NonLinearSolver', nls);

close all;
fnnobact = getPlotAfterStep(state0_comp_nobact, model_comp_nobact, schedule_comp, 'view', [0, 180], ...
    'field', 's:2');
[ws_comp_nobact, states_comp_nobact, report_nobact] = simulateScheduleAD(state0_comp_nobact, model_comp_nobact, schedule_comp, 'afterStepFn', fnnobact);

%% Run Compositional with Bacteria (if enabled)
fprintf('3. Running Compositional simulation (with bacteria)...\n');
close all;
fnbact = getPlotAfterStep(state0_comp_bact, model_comp_bact, schedule_comp, 'view', [0, 180], ...
    'field', 's:2');
[ws_comp_bact, states_comp_bact, report_bact] = simulateScheduleAD(state0_comp_bact, model_comp_bact, schedule_comp, 'afterStepFn', fnbact);

%% ============ PART 5: COMPARISON AND ANALYSIS =====================
fprintf('=== Comparing Results ===\n');

%% Extract Key Results
results_bo = extractResultsBO(states_bo, ws_bo, model_bo);
results_comp_nobact = extractResultsComp(states_comp_nobact, ws_comp_nobact, model_comp_nobact);
results_comp_bact = extractResultsComp(states_comp_bact, ws_comp_bact, model_comp_bact);

%% Calculate Performance Metrics
metrics_bo = calculateMetrics(results_bo, schedule_bo);
metrics_comp_nobact = calculateMetrics(results_comp_nobact, schedule_comp);
metrics_comp_bact = calculateMetrics(results_comp_bact, schedule_comp);

%% Plot Comparison
plotComparisonSummary(G, results_bo, results_comp_nobact, results_comp_bact, false);

%% Display Summary Table
fprintf('\n=== PERFORMANCE COMPARISON ===\n');
fprintf('Model                 | Avg Pressure | Max H2 Sat | Injected H2 | CPU Time\n');
fprintf('---------------------|--------------|------------|-------------|----------\n');
fprintf('Black-Oil            | %6.1f bar   | %6.4f    | %6.3e kg | %6.1f s\n', ...
    mean(results_bo.pressure)/barsa, max(results_bo.sG), ...
    sum(results_bo.H2_injected), report_bo.SimulationTime);

fprintf('Comp (No Bacteria)   | %6.1f bar   | %6.4f    | %6.3e kg | %6.1f s\n', ...
    mean(results_comp_nobact.pressure)/barsa, max(results_comp_nobact.sG), ...
    sum(results_comp_nobact.H2_injected), report_bo.SimulationTime);

fprintf('Comp (With Bacteria) | %6.1f bar   | %6.4f    | %6.3e kg | %6.1f s\n', ...
    mean(results_comp_bact.pressure)/barsa, max(results_comp_bact.sG), ...
    sum(results_comp_bact.H2_injected), report_bo.SimulationTime);

% Calculate bacterial effects
h2_loss = (max(results_comp_nobact.totMassH2) - max(results_comp_bact.totMassH2)) / ...
    max(results_comp_nobact.totMassH2) * 100;
fprintf('\nBacterial effects: H2 loss = %.2f%%\n', h2_loss);
fprintf('\n=== Simulation completed successfully ===\n');

%% ============ SUPPORTING FUNCTIONS =====================

function results = extractResultsBO(states, ws, model)
% Extract results from black-oil simulation
nT = numel(states);
results.pressure = zeros(nT, 1);
indH2 = 2;
results.yH2 = zeros(nT, 1);
results.sG = zeros(nT, 1);
results.H2_injected = zeros(nT, 1);
rhoOS = model.fluid.rhoOS;  % Density of oil-saturated phase
rhoGS = model.fluid.rhoGS;  % Density of gas-saturated phase
results.totMassH2 = zeros(nT, 1);

for i = 1:nT
    results.pressure(i) = mean(states{i}.pressure);
    results.sG(i) = max(states{i}.s(:,2));  % Gas saturation
    if i <= numel(ws) && ~isempty(ws{i})
        results.H2_injected(i) = max(0, ws{i}.qGs);  % Positive injection rates
    end
    results.yH2(i) =  max(rhoGS ./ (rhoGS + rhoOS .* states{i}.rv));
    if isfield(states{i}, 'FlowProps') && isfield(states{i}.FlowProps, 'ComponentTotalMass')
        results.totMassH2(i) = sum(states{i}.FlowProps.ComponentTotalMass{indH2});
    end
end

end

function results = extractResultsComp(states, ws, model)
% Extract results from compositional simulation
nT = numel(states);

% Get component indices
namecp = model.EOSModel.getComponentNames();
indH2 = find(strcmp(namecp, 'H2'));

results.pressure = zeros(nT, 1);
results.sG = zeros(nT, 1);
results.yH2 = zeros(nT, 1);
results.H2_injected = zeros(nT, 1);
results.totMassH2 = zeros(nT, 1);

for i = 1:nT
    results.pressure(i) = mean(states{i}.pressure);
    results.sG(i) = max(states{i}.s(:,2));
    results.yH2(i) = max(states{i}.y(:,indH2));

    if i <= numel(ws) && ~isempty(ws{i}) && isfield(ws{i}, 'H2')
        results.H2_injected(i) = max(0, ws{i}.H2);
    end

    % Total H2 mass in reservoir
    if isfield(states{i}, 'FlowProps') && isfield(states{i}.FlowProps, 'ComponentTotalMass')
        results.totMassH2(i) = sum(states{i}.FlowProps.ComponentTotalMass{indH2});
    end
end
end

function metrics = calculateMetrics(results, sched)
% Calculate performance metrics
metrics.avgPressure = mean(results.pressure);
metrics.maxPressure = max(results.pressure);
metrics.maxH2Saturation = max(results.sG);
metrics.totalInjectedH2 = sum(results.H2_injected);
end

function plotComparisonSummary(G, res_bo, res_comp_nobact, res_comp_bact, nobact)
% Plot comprehensive comparison
figure('Position', [100, 100, 1400, 900]);

% Pressure evolution
subplot(2,3,1);
plot(res_bo.pressure/barsa, 'b-', 'LineWidth', 2); hold on;
plot(res_comp_nobact.pressure/barsa, 'r-', 'LineWidth', 2);
if ~nobact
    plot(res_comp_bact.pressure/barsa, 'g-', 'LineWidth', 2);
    legend('Black-Oil', 'Comp (No Bact)', 'Comp (With Bact)', 'Location', 'best');
else
    legend('Black-Oil', 'Comp (No Bact)', 'Location', 'best');
end
xlabel('Time Step'); ylabel('Pressure (bar)');
title('Pressure Evolution'); grid on;

% H2 saturation
subplot(2,3,2);
plot(res_bo.sG, 'b-', 'LineWidth', 2); hold on;
plot(res_comp_nobact.sG, 'r-', 'LineWidth', 2);
if ~nobact
    plot(res_comp_bact.sG, 'g-', 'LineWidth', 2);
end
xlabel('Time Step'); ylabel('H2 Saturation');
title('Maximum H2 Saturation'); grid on;

% Cumulative H2 injection
subplot(2,3,3);
plot(cumsum(res_bo.H2_injected), 'b-', 'LineWidth', 2); hold on;
plot(cumsum(res_comp_nobact.H2_injected), 'r-', 'LineWidth', 2);
if ~nobact
    plot(cumsum(res_comp_bact.H2_injected), 'g-', 'LineWidth', 2);
end
xlabel('Time Step'); ylabel('Cumulative H2 (kg)');
title('Cumulative H2 Injection'); grid on;

% Final state comparison
subplot(2,3,4);
if ~nobact
    models = {'Black-Oil', 'Comp NoBact', 'Comp Bact'};
    pressures = [res_bo.pressure(end), res_comp_nobact.pressure(end), res_comp_bact.pressure(end)]/barsa;
    saturations = [res_bo.sG(end), res_comp_nobact.sG(end), res_comp_bact.sG(end)];
else
    models = {'Black-Oil', 'Comp NoBact'};
    pressures = [res_bo.pressure(end), res_comp_nobact.pressure(end)]/barsa;
    saturations = [res_bo.sG(end), res_comp_nobact.sG(end)];
end

yyaxis left;
bar(pressures);
ylabel('Final Pressure (bar)');
yyaxis right;
plot(saturations, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Final H2 Saturation');
set(gca, 'XTickLabel', models);
title('Final State Comparison'); grid on;

% H2 mole fraction
if isfield(res_comp_nobact, 'yH2')
    subplot(2,3,5);
    plot(res_bo.yH2, 'b-', 'LineWidth', 2); hold on;
    plot(res_comp_nobact.yH2, 'r-', 'LineWidth', 2); hold on;
    if ~nobact
        plot(res_comp_bact.yH2, 'g-', 'LineWidth', 2);
        legend('Black-Oil','Comp (No Bact)', 'Comp (With Bact)', 'Location', 'best');
    else
        legend('Comp (No Bact)', 'Location', 'best');
    end
    xlabel('Time Step'); ylabel('H2 Mole Fraction in Gas');
    title('H2 Gas Composition'); grid on;
end

% Total H2 mass in reservoir
if isfield(res_comp_nobact, 'totMassH2')
    subplot(2,3,6);

    plot(res_bo.totMassH2, 'b-', 'LineWidth', 2); hold on;
    plot(res_comp_nobact.totMassH2, 'r-', 'LineWidth', 2); hold on;
    if ~nobact
        plot(res_comp_bact.totMassH2, 'g-', 'LineWidth', 2);
        legend('Black-Oil','Comp (No Bact)', 'Comp (With Bact)', 'Location', 'best');
    end
    xlabel('Time Step'); ylabel('Total H2 Mass (kg)');
    title('H2 Mass in Reservoir'); grid on;
end

sgtitle('Hydrogen Storage: Black-Oil vs Compositional Model Comparison');
end