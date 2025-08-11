clear all
mrstModule add compositional ad-blackoil ad-core ad-props mrst-gui

% Set gravity and verbosity
gravity reset on
mrstVerbose = true;

%% Read the Eclipse deck file
deck = readEclipseDeck('C:\Users\elyesa\OneDrive - SINTEF\Documents\Projects\TUC_UHS_Benchmark\Simulation Cases/UHS_Benchmark_HighH2.DATA');

%% Prepare simulation parameters
bacteriamodel = true;
deck.PROPS.EOS = 'PR';
deck = convertDeckUnits(deck);

% Define compositional fluid
compFluid = TableCompositionalMixture({'water', 'Methane', 'Nitrogen', ...
    'CarbonDioxide', 'Ethane', 'Propane', 'n-butane', 'Hydrogen'}, ...
    {'H2O', 'C1', 'N2', 'CO2', 'C2', 'C3', 'NC4', 'H2'});

% Initialize Eclipse problem
[~, model, schedule, nls] = initEclipseProblemAD(deck, ...
    'getSchedule', true, 'getInitialState', false);

%% Set up equation of state
eos = SoreideWhitsonEquationOfStateModel([], compFluid, 'sw');
model.EOSModel = eos;
model.water = false;

%% Update fluid properties
[rhow, rhog] = deal(model.fluid.rhoWS, 1.481688 * kilogram/meter^3);
[viscow, viscog] = deal(model.fluid.muWr, 0.0094234 * centi * poise);
[cfw, cfg] = deal(model.fluid.cW, 8.1533e-3 / barsa);

fluid = initSimpleADIFluid('phases', 'OG', ...
    'mu', [viscow, viscog], ...
    'rho', [rhow, rhog], ...
    'pRef', 0 * barsa(), ...
    'c', [cfw, cfg], ...
    'n', [2, 2], ...
    'smin', [0.2, 0.1]);

% Assign relative permeability curves
deckBO = readEclipseDeck('\\wsl.localhost\Ubuntu\home\elyesa\Projects\MRST\modules\H2store\data\uhs_benchmark/UHS_BENCHMARK_RS.DATA');
fluid = assignSGOF(fluid, deckBO.PROPS.SGOF, struct('sat', 1, 'interp1d', @interpTable));
model.fluid = fluid;

%% Initialize state
T0 = deck.PROPS.TEMPVD{1}(2);
P0 = 81.6 * barsa();
s0 = [0.2 0.8];
z0 = [0.9023 0.07015 0.00405 0.020210 0.00272 0.0004 0.00015 0.00005];

% Set up backend
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);

if bacteriamodel
    compFluid = model.EOSModel.CompositionalMixture;
    arg = {model.G, model.rock, model.fluid, compFluid, ...
        false, diagonal_backend, 'eos', model.EOSModel, ...
        'oil', true, 'gas', true, ...
        'bacteriamodel', bacteriamodel, 'diffusioneffect', false, ...
        'liquidPhase', 'O', 'vaporPhase', 'G'};
    
    model = BiochemistryModel(arg{:});
    model.outputFluxes = false;
    nbact0 = 55;
    G = model.G;
    state0 = initCompositionalStateBacteria(model, P0.*ones(G.cells.num,1), ...
        T0, s0, z0, nbact0, eos);
else
    model.EOSModel = eos;
    state0 = initCompositionalState(model, P0.*ones(G.cells.num,1), T0, s0, z0);
end

% Adjust well compositions
for i = 1:length(schedule.control)
    schedule.control(i).W.compi = [0, 1];
end

%% Set up solvers
linsolve = nls.LinearSolver;
linsolve.amgcl_setup.block_size = 0;
linsolve.replaceInf = true;
linsolve.replaceNaN = true;
linsolve.decoupling = 'quasiIMPES';
%linsolve.maxIterations = 6;
nls = NonLinearSolver('LinearSolver', linsolve);
nls.maxIterations = 6;
nls.verbose = true;
model.verbose = true;

%% Run simulations
name = 'UHS_BENCHMARK_COMPOSITIONAL_BACT_TRUE';
problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);
simulatePackedProblem(problem);

% Run without bacteria for comparison
modelNoBact = model;
modelNoBact.bacteriamodel = false;
state0NoBact = state0;
state0NoBact.nbact = 0;

nameNoBact = 'UHS_BENCHMARK_COMPOSITIONAL_BACT_FALSE';
problemNoBact = packSimulationProblem(state0, modelNoBact, schedule, nameNoBact, 'NonLinearSolver', nls);
simulatePackedProblem(problemNoBact);

%% Post-process results
[~, states, reports] = getPackedSimulatorOutput(problem);
[~, statesNoBact, reportsNo] = getPackedSimulatorOutput(problemNoBact);

namecp = model.EOSModel.getComponentNames();
indH2 = find(strcmp(namecp, 'H2'));
indCO2 = find(strcmp(namecp, 'CO2'));
indCH4 = find(strcmp(namecp, 'C1'));

% Initialize arrays
nT = numel(states);
totalH2_bact = zeros(nT, 1);
totalH2_noBact = zeros(nT, 1);
totalCO2_bact = zeros(nT, 1);
totalCO2_noBact = zeros(nT, 1);
totalCH4_bact = zeros(nT, 1);
totalCH4_noBact = zeros(nT, 1);

for i = 1:nT
    totalH2_bact(i) = sum(states{i}.FlowProps.ComponentTotalMass{indH2});
    totalH2_noBact(i) = sum(statesNoBact{i}.FlowProps.ComponentTotalMass{indH2});
    totalCO2_bact(i) = sum(states{i}.FlowProps.ComponentTotalMass{indCO2});
    totalCO2_noBact(i) = sum(statesNoBact{i}.FlowProps.ComponentTotalMass{indCO2});
    totalCH4_bact(i) = sum(states{i}.FlowProps.ComponentTotalMass{indCH4});
    totalCH4_noBact(i) = sum(statesNoBact{i}.FlowProps.ComponentTotalMass{indCH4});
end

% Calculate percentages
H2_loss_percentage = ((totalH2_noBact - totalH2_bact) ./ totalH2_noBact) * 100;
CO2_loss_percentage = ((totalCO2_noBact - totalCO2_bact) ./ totalCO2_noBact) * 100;
CH4_loss_percentage = ((totalCH4_bact - totalCH4_noBact) ./ totalCH4_noBact) * 100;

% Display results
fprintf('Total H2 loss due to bacterial effects: %.2f%%\n', H2_loss_percentage(end));