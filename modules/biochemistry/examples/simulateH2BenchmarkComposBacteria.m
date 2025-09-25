%% 3D Benchmark of Hydrogen Storage Cycling in an Aquifer
% ===========================================================================
% This benchmark simulates cyclic hydrogen storage in an aquifer using 
% compositional models with bacterial effects.
%
% Original Benchmark:
% Hogeweg, S., Strobel, G., & Hagemann, B. (2022). Benchmark study for the 
% simulation of underground hydrogen storage operations. Comput Geosci, 26, 1367–1378.
% https://doi.org/10.1007/s10596-022-10163-5
%
% Reference for the case using blackoil model: Ahmed, E., et al., 2024, "Phase behavior and black-oil simulations
% of hydrogen storage in saline aquifers."
% Reference for the current case with fully compositional and biochmestry model: https://www.sciencedirect.com/science/article/pii/S0360319925039473
%
% Key aspcets:
% - Initial brine-filled reservoir
% - SW EoS is used
% - Brooks–Corey capillary pressure model
% - Compositional model with 8 components
% - Comparison of bacterial vs non-bacterial effects
% ---------------------------------------------------------------------------

clear; clc;
mrstModule add compositional ad-blackoil ad-core ad-props deckformat mrst-gui h2store biochemistry

gravity reset on
mrstVerbose = true;
%% Notice on Computational Cost
warning('ComputationalCost:High', ...
       ['Please be advised that this example takes a long time ', ...
        'to run']);

%% ============ DOWNLOAD AND SETUP BENCHMARK DATA =====================
fprintf('=== Setting up Benchmark Case ===\n');

benchmarkZipUrl = 'https://www.ite.tu-clausthal.de/fileadmin/ITE/images/Research/Energy_and_Gas_Storage/TUC_UHS_Benchmark_v1.zip';
zipFile = 'TUC_UHS_Benchmark_v1.zip';

% Download zip archive
if ~isfile(zipFile)
    fprintf('Downloading benchmark data...\n');
    websave(zipFile, benchmarkZipUrl);
end

% Extract to folder
if ~isfolder('TUC_UHS_Benchmark_v1')
    unzip(zipFile, 'TUC_UHS_Benchmark_v1');
end

% Find .DATA file
files = dir(fullfile('TUC_UHS_Benchmark_v1', '**', '*.DATA'));
if ~isempty(files)
    mainDataFile = fullfile(files(1).folder, files(1).name);
    fprintf('Benchmark ready: %s\n', mainDataFile);
else
    error('No .DATA file found in extracted archive.');
end

% Black-oil data for rel-perm
blackoilDataFile = fullfile(getDatasetPath('uhs_benchmark'), 'UHS_BENCHMARK_RS.DATA');

%% ============ LOAD ECLIPSE DECKS =====================
deck = readEclipseDeck(mainDataFile);
deck.PROPS.EOS = 'PR';
deck = convertDeckUnits(deck);

deckBO = readEclipseDeck(blackoilDataFile);
deckBO = convertDeckUnits(deckBO);

%% ============ COMPOSITIONAL FLUID DEFINITION =====================
compFluid = TableCompositionalMixture({'water', 'Methane', 'Nitrogen', ...
    'CarbonDioxide', 'Ethane', 'Propane', 'n-butane', 'Hydrogen'}, ...
    {'H2O', 'C1', 'N2', 'CO2', 'C2', 'C3', 'NC4', 'H2'});

%% ============ MODEL INITIALIZATION =====================
[~, model, schedule, nls_bact] = initEclipseProblemAD(deck, 'getSchedule', true, 'getInitialState', false);

% Equation of state
eos = EquationOfStateModel([], compFluid, 'sw');
model.EOSModel = eos;
model.water = false;

%% ============ FLUID PROPERTIES SETUP =====================
[rhow, rhog] = deal(model.fluid.rhoWS, 1.481688 * kilogram/meter^3);
[viscow, viscog] = deal(model.fluid.muWr, 0.0094234 * centi * poise);
[cfw, cfg] = deal(model.fluid.cW, 8.1533e-3 / barsa);

fluid = initSimpleADIFluid('phases', 'OG', ...
    'mu', [viscow, viscog], 'rho', [rhow, rhog], 'pRef', 0*barsa(), ...
    'c', [cfw, cfg], 'n', [2, 2], 'smin', [0.2, 0.1]);

% Assign relative permeability
fluid = assignSGOF(fluid, deckBO.PROPS.SGOF, struct('sat', 1, 'interp1d', @interpTable));
model.fluid = fluid;

%% ============ INITIAL CONDITIONS =====================
T0 = deck.PROPS.TEMPVD{1}(2);
P0 = 81.6 * barsa();
s0 = [0.2, 0.8];
z0 = [0.9023, 0.07015, 0.00405, 0.020210, 0.00272, 0.0004, 0.00015, 0.00005];

backend = DiagonalAutoDiffBackend('modifyOperators', true);
G = model.G;

%% ============ SETUP MODELS WITH SEPARATE SOLVERS =====================
fprintf('\n=== Setting up Models ===\n');

% Common model arguments
arg_common = {G, model.rock, model.fluid, compFluid, false, backend, ...
    'water', false, 'oil', true, 'gas', true, ...
    'liquidPhase', 'O', 'vaporPhase', 'G'};

% Model 1: With Bacteria
fprintf('1. Setting up model WITH bacteria...\n');
model_bact = BiochemistryModel(arg_common{:}, 'bacteriamodel', true);
model_bact.outputFluxes = false;

% Initialize state with bacteria
nbact0 = 55;
state0_bact = initCompositionalStateBacteria(model_bact, P0*ones(G.cells.num,1), T0, s0, z0, nbact0, eos);

% Solver for bacterial model
lsolve_bact = nls_bact.LinearSolver;% selectLinearSolverAD(model_bact);
lsolve_bact.amgcl_setup.block_size = 0;
lsolve_bact.decoupling = 'quasiIMPES';
nls_bact = NonLinearSolver('LinearSolver', lsolve_bact);
nls_bact.maxIterations = 12;
nls_bact.timeStepSelector.maxTimestep = 5*day;
model_bact.verbose = true;

% Model 2: Without Bacteria  
fprintf('2. Setting up model WITHOUT bacteria...\n');
model_nobact = BiochemistryModel(arg_common{:}, 'bacteriamodel', false);

% Initialize state without bacteria
state0_nobact = initCompositionalState(model_nobact, P0*ones(G.cells.num,1), T0, s0, z0);

% Solver for non-bacterial model
lsolve_nobact = selectLinearSolverAD(model_nobact);
lsolve_nobact.decoupling = 'quasiIMPES';
nls_nobact = NonLinearSolver('LinearSolver', lsolve_nobact);
model_nobact.verbose = true;

% Adjust well compositions
for i = 1:length(schedule.control)
    schedule.control(i).W.compi = [0, 1];
end

%% ============ RUN SIMULATIONS =====================
fprintf('\n=== Running Simulations ===\n');

% Run with bacteria first
fprintf('1. Running WITH bacteria...\n');
problem_bact = packSimulationProblem(state0_bact, model_bact, schedule, ...
    'UHS_BENCHMARK_WITH_BACTERIA', 'NonLinearSolver', nls_bact);
simulatePackedProblem(problem_bact);
fprintf('   ✓ Completed with bacteria\n');

% Run without bacteria
fprintf('2. Running WITHOUT bacteria...\n');
problem_nobact = packSimulationProblem(state0_nobact, model_nobact, schedule, ...
    'UHS_BENCHMARK_NO_BACTERIA', 'NonLinearSolver', nls_nobact);
simulatePackedProblem(problem_nobact);
fprintf('   ✓ Completed without bacteria\n');

%% ============ POST-PROCESSING =====================
fprintf('\n=== Post-processing Results ===\n');

% Load results
[ws_bact, states_bact, ~] = getPackedSimulatorOutput(problem_bact);
[ws_nobact, states_nobact, ~] = getPackedSimulatorOutput(problem_nobact);

% Get component indices
namecp = model_bact.EOSModel.getComponentNames();
indH2  = find(strcmp(namecp, 'H2'));
indCO2 = find(strcmp(namecp, 'CO2'));
indCH4 = find(strcmp(namecp, 'C1'));

nT = numel(states_bact);

% Initialize result arrays
totalH2_bact   = zeros(nT, 1); totalH2_nobact   = zeros(nT, 1);
totalCO2_bact  = zeros(nT, 1); totalCO2_nobact  = zeros(nT, 1);
totalCH4_bact  = zeros(nT, 1); totalCH4_nobact  = zeros(nT, 1);

% Extract component masses
for i = 1:nT
    totalH2_bact(i)   = sum(states_bact{i}.FlowProps.ComponentTotalMass{indH2});
    totalH2_nobact(i) = sum(states_nobact{i}.FlowProps.ComponentTotalMass{indH2});
    totalCO2_bact(i)  = sum(states_bact{i}.FlowProps.ComponentTotalMass{indCO2});
    totalCO2_nobact(i)= sum(states_nobact{i}.FlowProps.ComponentTotalMass{indCO2});
    totalCH4_bact(i)  = sum(states_bact{i}.FlowProps.ComponentTotalMass{indCH4});
    totalCH4_nobact(i)= sum(states_nobact{i}.FlowProps.ComponentTotalMass{indCH4});
end

% Calculate effects
H2_loss  = ((totalH2_nobact - totalH2_bact) ./ totalH2_nobact) * 100;
CO2_loss = ((totalCO2_nobact - totalCO2_bact) ./ totalCO2_nobact) * 100;
CH4_gain = ((totalCH4_bact - totalCH4_nobact) ./ totalCH4_nobact) * 100;

%% ============ RESULTS SUMMARY =====================
fprintf('\n=== BACTERIAL EFFECTS SUMMARY ===\n');
fprintf('H2 loss due to bacteria:  %.2f%%\n', H2_loss(end));
fprintf('CO2 loss due to bacteria: %.2f%%\n', CO2_loss(end));
fprintf('CH4 production due to bacteria: %.2f%%\n', CH4_gain(end));

fprintf('\n=== Final Mass Balances ===\n');
fprintf('               WITH Bacteria  WITHOUT Bacteria\n');
fprintf('H2 mass:       %12.3e  %12.3e kg\n', totalH2_bact(end), totalH2_nobact(end));
fprintf('CO2 mass:      %12.3e  %12.3e kg\n', totalCO2_bact(end), totalCO2_nobact(end));
fprintf('CH4 mass:      %12.3e  %12.3e kg\n', totalCH4_bact(end), totalCH4_nobact(end));

fprintf('\n=== Benchmark simulation completed ===\n');

%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>