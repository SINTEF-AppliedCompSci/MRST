%% Test: Molecular diffusion ON/OFF for 3 scenarios (clog, noclog, nobact)
% Runs 6 packed simulations and then loads + plots results.
%
% This script compares the effect of molecular diffusion on hydrogen storage
% simulations with three different scenarios:
%   A: Bacteria + clogging (dynamic porosity/permeability)
%   B: Bacteria only (no clogging, static rock properties)
%   C: No bacteria (baseline case)
%
% For each scenario, simulations are run with molecular diffusion ON and OFF.
%
% SEE ALSO:%   `BiochemistryModel`, `setupBioCloggingModel`, `runCase`
%
%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

MRST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST. If not, see <http://www.gnu.org/licenses/>.
%}

clearvars;

mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui upr test-suite spe10
mrstModule add compositional h2-biochem h2store

baseName = 'H2_STORAGE_DOME_TRAP';
dataPath = getDatasetPath('h2storage');
dataFile = fullfile(dataPath, 'H2STORAGE_RS.DATA');
deck = readEclipseDeck(dataFile);

warning('ComputationalCost:High', 'Running multiple cases; reduce cycles if needed.');

%% Base black-oil setup
[~, ~, state0Bo, modelBo, scheduleBo, ~] = modelForSimple2DAquifer(deck,'numcycles',10);

%% Convert BO -> compositional (capture G/rock/fluid explicitly)
tmp = convertBlackOilModelToCompositionalModel(modelBo);
G = tmp.G;
rock0  = tmp.rock;
fluid0 = tmp.fluid;

state0Base = convertBlackOilStateToCompositional(modelBo, state0Bo);

%% Mixture + EOS
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, ...
    {'H2O',   'H2',      'CO2',           'C1'});
EOS = EquationOfStateModel([], compFluid, 'sw');

nc = G.cells.num;
T0 = 273.15 + 44.35;
comp0 = repmat([0.8480, 1.0e-5, 1.0e-5, 0.1530], nc, 1);

%% Schedule edits (reuse for all cases)
schedule = scheduleBo;
for i = 1:numel(schedule.control)
    schedule.control(i).W.compi = [0, 1];
    if strcmp(schedule.control(i).W.name, 'cushion') && i < 11
        schedule.control(i).W.components = [0.0, 0.1, 0.9, 0.0];
    else
        schedule.control(i).W.components = [0.0, 0.95, 0.05, 0.0];
    end
    schedule.control(i).W.T = T0;
    schedule.control(i).bc = [];
end

%% Nonlinear solver (shared)
nls = NonLinearSolver();

%% Scenario A: bacteria + clogging (dynamic poro/perm)
nbact0  = 5;
nc_bact = 120;
cp      = 1.0;

tmp2 = tmp;
[tmp2, poroA, permA] = setupBioCloggingModel(tmp2, nbact0, nc_bact, cp);
rockA  = tmp2.rock;
fluidA = tmp2.fluid;

%% Scenario B/C: static rock/fluid baseline (no clogging)
rockStatic = rock0;
rockStatic.poro = poroA;
rockStatic.perm = permA;

fluidStatic = fluid0;
fluidStatic.pvMultR = @(p, nbact) 1;

%% Run all 6 cases
pA_off = runCase([baseName '_A_CLOG_MD_OFF'],  G, rockA,      fluidA,      compFluid, EOS, modelBo.gravity, state0Base, comp0, schedule, nls, true,  false, nbact0);
pA_on  = runCase([baseName '_A_CLOG_MD_ON'],   G, rockA,      fluidA,      compFluid, EOS, modelBo.gravity, state0Base, comp0, schedule, nls, true,  true,  nbact0);

pB_off = runCase([baseName '_B_NOCLG_MD_OFF'], G, rockStatic, fluidStatic, compFluid, EOS, modelBo.gravity, state0Base, comp0, schedule, nls, true,  false, nbact0);
pB_on  = runCase([baseName '_B_NOCLG_MD_ON'],  G, rockStatic, fluidStatic, compFluid, EOS, modelBo.gravity, state0Base, comp0, schedule, nls, true,  true,  nbact0);

%pC_off = runCase([baseName '_C_NOBAC_MD_OFF'], G, rockStatic, fluidStatic, compFluid, EOS, modelBo.gravity, state0Base, comp0, schedule, nls, false, false, 0);
pC_on  = runCase([baseName '_C_NOBAC_MD_ON'],  G, rockStatic, fluidStatic, compFluid, EOS, modelBo.gravity, state0Base, comp0, schedule, nls, false, true,  0);

%% Extract results
[wsA_off, stA_off] = getPackedSimulatorOutput(pA_off);
[wsA_on,  stA_on ] = getPackedSimulatorOutput(pA_on);

[wsB_off, stB_off] = getPackedSimulatorOutput(pB_off);
[wsB_on,  stB_on ] = getPackedSimulatorOutput(pB_on);

[wsC_off, stC_off] = getPackedSimulatorOutput(pC_off);
[wsC_on,  stC_on ] = getPackedSimulatorOutput(pC_on);

%% Plot (simple)
figure('Name','A (clog) MD OFF','Color','w'); plotToolbar(pA_off.SimulatorSetup.model.G, stA_off); title('A (clog) MD OFF','Interpreter','none');
figure('Name','A (clog) MD ON','Color','w');  plotToolbar(pA_on.SimulatorSetup.model.G,  stA_on ); title('A (clog) MD ON','Interpreter','none');

figure('Name','B (no clog) MD OFF','Color','w'); plotToolbar(pB_off.SimulatorSetup.model.G, stB_off); title('B (no clog) MD OFF','Interpreter','none');
figure('Name','B (no clog) MD ON','Color','w');  plotToolbar(pB_on.SimulatorSetup.model.G,  stB_on ); title('B (no clog) MD ON','Interpreter','none');

figure('Name','C (no bact) MD OFF','Color','w'); plotToolbar(pC_off.SimulatorSetup.model.G, stC_off); title('C (no bact) MD OFF','Interpreter','none');
figure('Name','C (no bact) MD ON','Color','w');  plotToolbar(pC_on.SimulatorSetup.model.G,  stC_on ); title('C (no bact) MD ON','Interpreter','none');

figure('Name','Wells A OFF vs ON','Color','w'); plotWellSols({wsA_off, wsA_on}); legend({'OFF','ON'},'Location','best'); title('Wells A (clog) OFF vs ON','Interpreter','none');
figure('Name','Wells B OFF vs ON','Color','w'); plotWellSols({wsB_off, wsB_on}); legend({'OFF','ON'},'Location','best'); title('Wells B (no clog) OFF vs ON','Interpreter','none');
figure('Name','Wells C OFF vs ON','Color','w'); plotWellSols({wsC_off, wsC_on}); legend({'OFF','ON'},'Location','best'); title('Wells C (no bact) OFF vs ON','Interpreter','none');

%% Optional: save problems so you can plot later without rerunning
save('mdiffPackedProblems.mat', 'pA_off','pA_on','pB_off','pB_on','pC_off','pC_on');

function problem = runCase(caseName, G, rock, fluid, compFluid, EOS, gravity, state0Base, comp0, schedule, nls, bacteriamodel, mdiff, nbact0)
fprintf('Running case: %s\n', caseName);

backend = DiagonalAutoDiffBackend('modifyOperators', true);

model = BiochemistryModel( ...
    G, rock, fluid, compFluid, ...
    false, backend, ...
    'oil', true, 'gas', true, ...
    'bacteriamodel', bacteriamodel, ...
    'molecularDiffusion', mdiff, ...
    'liquidPhase', 'O', 'vaporPhase', 'G');

model.EOSModel = EOS;
model.gravity  = gravity;
model.outputFluxes = false;

% Diagnostic: verify diffusion is actually wired into ComponentTotalFlux
try
    ctf = model.FlowDiscretization.getStateFunction('ComponentTotalFlux');
    fprintf('  [diag] bacteriamodel=%d mdiff=%d: FlowDisc=%s, ComponentTotalFlux=%s\n', ...
        bacteriamodel, mdiff, class(model.FlowDiscretization), class(ctf));
catch me
    fprintf('  [diag] Could not query flux SF: %s\n', me.message);
end

% Initial state
T0 = schedule.control(1).W.T;
if bacteriamodel
    state0 = initCompositionalStateBacteria(model, state0Base.pressure, T0, state0Base.s, comp0, nbact0, EOS);
else
    state0 = initCompositionalState(model, state0Base.pressure, T0, state0Base.s, comp0, EOS);
    if isfield(state0, 'nbact')
        state0.nbact(:) = 0;
    end
end

% Linear solver
lsolve = selectLinearSolverAD(model);
nls.LinearSolver = lsolve;

% Pack + run
problem = packSimulationProblem(state0, model, schedule, caseName, 'NonLinearSolver', nls);
simulatePackedProblem(problem, 'RestartStep', 1);
end