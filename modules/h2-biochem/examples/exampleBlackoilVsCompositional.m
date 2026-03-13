%% Comprehensive Hydrogen Storage: Black-Oil vs Compositional Comparison
% ===========================================================================
% This script compares three modeling approaches:
% 1. Black-oil model with tabulated PVT data (from h2store module)
% 2. Compositional model without bacterial effects
% 3. Compositional model with bacterial growth
%
% SEE ALSO:
%   `simpleExampleBlackOilPVT`, `simpleBlackOil2DModel`
%
%
%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

%clearvars; clc;
mrstModule add ad-blackoil ad-core ad-props mrst-gui deckformat
mrstModule add compositional h2-biochem h2store

%% ============ PART 1: BLACK-OIL MODEL (Original Script) =====================
fprintf('=== Setting up Black-Oil Model ===\n');
simpleBlackOil2DModel();

%% ============ PART 2: COMPOSITIONAL MODELS =====================
fprintf('=== Setting up Compositional Models ===\n');

%% Define Compositional Fluid (H2, H2O, CO2, CH4)
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, ...
    {'H2O', 'H2', 'CO2', 'C1'});

%% Common Parameters for Compositional Models
z0 = [1.0, 0.0, 0.00, 0.0];  % Initial composition [H2O, H2, CO2, C1]
T0 =  273.15 +60;
%% Model 1: Compositional without Bacteria
backend = DiagonalAutoDiffBackend('modifyOperators', true);
model_comp_nobact = BiochemistryModel(G, rock, fluid_bo, compFluid, ...
    true, backend, 'water', false, 'oil', true, 'gas', true, ...
    'bacteriamodel', false, 'liquidPhase', 'O', 'vaporPhase', 'G');
model_comp_nobact.gravity = grav;
model_comp_nobact.nonlinearTolerance = 1.0e-6;
model_comp_nobact.OutputStateFunctions{end+1} = 'ComponentPhaseMass';
% Initialize compositional state without bacteria
state0_comp_nobact = initCompositionalState(model_comp_nobact, pInit, T0, [1.0, 0.0], z0);

%% Model 2: Compositional with Bacteria (if enabled)
model_comp_bact = BiochemistryModel(G, rock, fluid_bo, compFluid, ...
    true, backend, 'water', false, 'oil', true, 'gas', true, ...
    'bacteriamodel', true, 'liquidPhase', 'O', 'vaporPhase', 'G');
model_comp_bact.gravity = grav;
model_comp_bact.OutputStateFunctions{end+1} = 'ComponentPhaseMass';
model_comp_bact.nonlinearTolerance = 1.0e-6;
nbact0 = 100;
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

%% Run Black-Oil Simulation
fprintf('1. Running Black-Oil simulation...\n');
clf;
fn = getPlotAfterStep(state0_bo, model_bo, schedule_bo, ...
    'view', [0, 180], 'field', 's:2');
[ws_bo, states_bo, report_bo] = simulateScheduleAD(state0_bo, model_bo, schedule_bo, 'afterStepFn', fn);

%% Run Compositional without Bacteria
fprintf('2. Running Compositional simulation (no bacteria)...\n');

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
results_bo = extractResults(states_bo, ws_bo, model_bo);
results_comp_nobact = extractResults(states_comp_nobact, ws_comp_nobact, model_comp_nobact);
results_comp_bact = extractResults(states_comp_bact, ws_comp_bact, model_comp_bact);

%% Calculate Performance Metrics
metrics_bo = calculateMetrics(results_bo);
metrics_comp_nobact = calculateMetrics(results_comp_nobact);
metrics_comp_bact = calculateMetrics(results_comp_bact);

%% Plot Comparison
plotComparisonSummary(results_bo, results_comp_nobact, results_comp_bact);

% Calculate H2 loss due to bacteria (comparing compositional with vs without bacteria)
% Reference: compositional without bacteria
h2_loss_due_to_bacteria = (max(results_comp_nobact.totMassH2) - max(results_comp_bact.totMassH2)) / ...
    max(results_comp_nobact.totMassH2) * 100;
fprintf('\nH2 loss due to bacteria: %.2f%%\n', h2_loss_due_to_bacteria);

% Calculate H2 mass error between blackoil and compositional without bacteria
% Reference: compositional without bacteria
h2_error_blackoil_vs_comp_nobact = (max(results_bo.totMassH2) - max(results_comp_nobact.totMassH2)) / ...
    max(results_comp_nobact.totMassH2) * 100;
fprintf('\nH2 mass error - Blackoil vs Compositional (no bacteria): %.2f%%\n', h2_error_blackoil_vs_comp_nobact);

% Calculate H2 mass error between blackoil and compositional with bacteria
% Reference: compositional with bacteria
h2_error_blackoil_vs_comp_bact = (max(results_bo.totMassH2) - max(results_comp_bact.totMassH2)) / ...
    max(results_comp_bact.totMassH2) * 100;
fprintf('\nH2 mass error - Blackoil vs Compositional (with bacteria): %.2f%%\n', h2_error_blackoil_vs_comp_bact);

fprintf('\n=== Simulation completed successfully ===\n');


%% ============ SUPPORTING FUNCTIONS =====================

function plotComparisonSummary(res_bo, res_comp_nobact, res_comp_bact)
% Plot comparison between Black-Oil, Compositional (No Bacteria), and Compositional (With Bacteria)

figure('Position', [100, 100, 1200, 800]);

res_list = {res_bo, res_comp_nobact, res_comp_bact};
labels   = {'Black-Oil', 'Comp (No Bact)', 'Comp (With Bact)'};
colors   = {'b-', 'r-', 'g-'};

% Create subplots
subplot(2,3,1);
plotSeries(res_list, 'pressure', 'Pressure (bar)', 'Pressure Evolution', false, colors);
subplot(2,3,2);
plotSeries(res_list, 'sG', 'H2 Saturation', 'Mean H2 Saturation', false, colors);
subplot(2,3,3);
plotSeries(res_list, 'H2_diss', ' H2 Mass in Water(kg)', 'Dissolved H2', true, colors);

% Final state comparison
subplot(2,3,4);
pressures   = cellfun(@(r) r.pressure(end)/barsa, res_list);
saturations = cellfun(@(r) r.sG(end), res_list);
yyaxis left; bar(pressures); ylabel('Final Pressure (bar)');
yyaxis right; plot(saturations, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
set(gca, 'XTickLabel', labels); title('Final State Comparison'); grid on;

% H2 mole fraction
if isfield(res_comp_nobact, 'yH2')
    subplot(2,3,5);
    plotSeries(res_list, 'yH2', 'H2 Mole Fraction in Gas', 'Mean yH2', false, colors);
end

% Total H2 mass
if isfield(res_comp_nobact, 'totMassH2')
    subplot(2,3,6);
    plotSeries(res_list, 'totMassH2', 'Total H2 Mass (kg)', 'H2 Mass in Reservoir', false, colors);
end

legend(labels, 'Location', 'bestoutside');
sgtitle('Hydrogen Storage: Black-Oil vs Compositional Model Comparison');

end

function plotSeries(res_list, field, ylab, ttl, cumulative, colors)
% Helper function to plot series for multiple results
hold on; grid on;
for i = 1:numel(res_list)
    if isfield(res_list{i}, field)
        y = res_list{i}.(field);
        if cumulative
            y = cumsum(y);
        end
        plot(y, colors{i}, 'LineWidth', 2);
    end
end
xlabel('Time Step'); ylabel(ylab); title(ttl);
end

function results = extractResults(states, ws, model)
% Extract results from simulation (works for both black-oil and compositional models)
%
% Parameters:
%   states - Cell array of state variables
%   ws     - Well solution data
%   model  - Simulation model (GenericBlackOilModel or BiochemistryModel)

nT = numel(states);
results.pressure     = zeros(nT, 1);
results.yH2          = zeros(nT, 1);
results.sG           = zeros(nT, 1);
results.H2_diss  = zeros(nT, 1);
results.totMassH2    = zeros(nT, 1);

% Determine model type and setup indices
isBlackOil      = isa(model, 'GenericBlackOilModel');
isCompositional = isa(model, 'BiochemistryModel');

if isBlackOil
    % Black-oil model setup
    rhoOS = model.fluid.rhoOS;
    rhoGS = model.fluid.rhoGS;
    names = model.getPhaseNames();
    nph = model.getNumberOfPhases;
    for ph = 1:nph
        switch names(ph)
            case 'O'
                L_ix = ph;
            case 'G'
                V_ix = ph;
        end
    end
    compIndex = V_ix; % Gas component index for total mass
elseif isCompositional
    % Compositional model setup
    namecp   = model.EOSModel.getComponentNames();
    indH2    = find(strcmp(namecp, 'H2'));
    V_ix     = model.getVaporIndex();
    L_ix     = model.getLiquidIndex();
    compIndex = indH2; % H2 index in compositional model
else
    error('Unknown model type: %s', class(model));
end

% Loop over timesteps
for i = 1:nT
    state = states{i};

    % Common quantities
    results.pressure(i) = mean(state.pressure);
    results.sG(i)       = mean(state.s(:, V_ix));

    if isBlackOil
        results.yH2(i) = mean(rhoGS .*state.rs./ (rhoGS + rhoOS .* state.rv));
        results.H2_diss(i) =sum(state.FlowProps.ComponentPhaseMass{compIndex,L_ix});
    else
        results.yH2(i) = mean(state.y(:, indH2));
        results.H2_diss(i) = sum(state.FlowProps.ComponentPhaseMass{compIndex,L_ix});

    end

    % Total H2 mass (if available)
    if isfield(state, 'FlowProps') && isfield(state.FlowProps, 'ComponentTotalMass')
        if length(state.FlowProps.ComponentTotalMass) >= compIndex
            results.totMassH2(i) = sum(state.FlowProps.ComponentTotalMass{compIndex});
        end
    end
end
end

function metrics = calculateMetrics(results)
% Calculate performance metrics
metrics.avgPressure = mean(results.pressure);
metrics.maxPressure = max(results.pressure);
metrics.maxH2Saturation = max(results.sG);
metrics.H2diss = sum(results.H2_diss);
if isfield(results, 'totMassH2')
    metrics.finalH2Mass = results.totMassH2(end);
end
end