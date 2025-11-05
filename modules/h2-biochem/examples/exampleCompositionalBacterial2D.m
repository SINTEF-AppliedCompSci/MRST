%% 2D Compositional Hydrogen Storage with Bio-Clogging
% ===========================================================================
% This MRST example simulates the storage of hydrogen in a saline aquifer
% and investigates the effects of bacterial on reservoir performance.
%
% Key features:
%  - 2D dome-shaped aquifer with caprock and bedrock layers.
%  - Injection of hydrogen into pure water and brine with salinity.
%  - Bio-clogging effects modeled via porosity and permeability reduction.
%  - Compositional model with H2, H2O, CO2, and CH4 components.
%  - Scenarios:
%       1) With bacterial effects and clogging
%       2) With bacterial effects but without clogging
%       3) Abiotic (no bacteria)
%  - EOS: Soreide-Whitson
%
% Reference for original case: Ahmed, E., et al., 2024, "Phase behavior and black-oil simulations
% of hydrogen storage in saline aquifers."
% Reference for the current case: https://www.sciencedirect.com/science/article/pii/S0360319925039473
% ---------------------------------------------------------------------------

clearvars;
mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui upr test-suite spe10
mrstModule add compositional h2-biochem h2store
%% Define case identifiers and Eclipse deck
baseName = 'H2_STORAGE_DOME_TRAP';
dataPath = getDatasetPath('h2storage');
dataFile = fullfile(dataPath, 'H2STORAGE_RS.DATA');
deck = readEclipseDeck(dataFile);

%% Warn about computational cost
warning('ComputationalCost:High', ...
    ['This is a 10-cycle injection example; consider reducing cycles for faster runs.']);

%% Set up black-oil model and schedule
[~, ~, state0Bo, modelBo, scheduleBo, ~] = modelForSimple2DAquifer(deck);

%% Convert black-oil to compositional model
model = convertBlackOilModelToCompositionalModel(modelBo);
state0 = convertBlackOilStateToCompositional(modelBo, state0Bo);

%% Define compositional fluid (H2, H2O, CO2, CH4)
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, ...
    {'H2O', 'H2', 'CO2', 'C1'});
EOS = EquationOfStateModel([], compFluid, 'sw');
model.EOSModel = EOS;

nc = model.G.cells.num;
T0 = 273.15 + 44.35;  % Initial temperature in K
comp0 = repmat([0.8480, 1.0e-5, 1.0e-5, 0.1530], nc, 1);

%% Bio-clogging parameters
bacteriamodel = true;
clogModel = true;
if clogModel
    nbact0 = 5; % Initial bacteria
    nc_bact = 120;
    cp = 1.0;   % Clogging coefficient
    [model, poro0, perm0] = setupBioCloggingModel(model, nbact0, nc_bact, cp);
    caseNameWithClogging = [baseName '_WITH_CLOGGING'];
else
    poro0 = model.rock.poro;
    perm0 = model.rock.perm;
end
%% Setup model
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
arg = {model.G, model.rock, model.fluid, compFluid,...
    false, diagonal_backend, 'oil', true, 'gas', true, ... % Define phases for water-oil system
    'bacteriamodel', bacteriamodel, 'liquidPhase', 'O', ...
    'vaporPhase', 'G'}; % Set phases and EOS model
model = BiochemistryModel(arg{:});
model.gravity = modelBo.gravity;

%% Initialize compositional state
if bacteriamodel
    state0 = initCompositionalStateBacteria(model, state0.pressure, T0, state0.s, comp0, nbact0, EOS);
else
    state0 = initCompositionalState(model, state0.pressure, T0, state0.s, state0.components, EOS);
end

%% Update schedule controls for compositional injection
schedule = scheduleBo;
for i = 1:numel(schedule.control)
    schedule.control(i).W.compi = [0, 1]; % Well components
    if strcmp(schedule.control(i).W.name, 'cushion') && i<11
        schedule.control(i).W.components = [0.0, 0.1, 0.9, 0.0];
    else
        schedule.control(i).W.components = [0.0, 0.95, 0.05, 0.0];
    end
    schedule.control(i).W.T = T0; % Temperature
    schedule.control(i).bc = [];   % Remove BCs from controls
end
model.outputFluxes = false;

%% Plot initial grid, porosity, and permeability
figure;
subplot(1,2,1);
plotCellData(model.G, poro0); hold on;
plotGrid(model.G, schedule.control(1).W.cells, 'FaceColor','red','LineStyle','none');
title('Porosity'); axis off tight;

subplot(1,2,2);
plotCellData(model.G, log10(perm0(:,1))); hold on;
plotGrid(model.G, schedule.control(1).W.cells, 'FaceColor','red','LineStyle','none');
title('Permeability (log10)'); axis off tight;
sgtitle('2D Dome-Shaped Aquifer');

%% Initialize solvers
nls = NonLinearSolver();
lsolve = selectLinearSolverAD(model);
nls.LinearSolver = lsolve;

%% Pack and run simulations
% --- With bio-clogging
problemWithClogging = packSimulationProblem(state0, model, schedule, caseNameWithClogging, 'NonLinearSolver', nls);
simulatePackedProblem(problemWithClogging);

% --- Without clogging
modelNoClogging = model;
modelNoClogging.rock.perm = perm0;
modelNoClogging.rock.poro = poro0;
modelNoClogging.fluid.pvMultR = @(p, nbact) 1;
modelNoClogging = BiochemistryModel(modelNoClogging.G, modelNoClogging.rock, modelNoClogging.fluid, compFluid, ...
    false, DiagonalAutoDiffBackend('modifyOperators', true), 'oil', true, 'gas', true, ...
    'bacteriamodel', true, 'liquidPhase', 'O', 'vaporPhase', 'G');
state0NoClogging = state0;
caseNameNoClogging = [baseName '_NO_CLOGGING'];
problemNoClogging = packSimulationProblem(state0NoClogging, modelNoClogging, schedule, caseNameNoClogging, 'NonLinearSolver', nls);
simulatePackedProblem(problemNoClogging);

% --- Without bacteria
modelNoBact = model;
modelNoBact.rock.perm = perm0;
modelNoBact.rock.poro = poro0;
modelNoBact.fluid.pvMultR = @(p, nbact) 1;
modelNoBact = BiochemistryModel(modelNoBact.G, modelNoBact.rock, modelNoBact.fluid, compFluid, ...
    false, DiagonalAutoDiffBackend('modifyOperators', true), 'oil', true, 'gas', true, ...
    'bacteriamodel', false, 'liquidPhase', 'O', 'vaporPhase', 'G');
state0NoBact = state0;
state0NoBact.nbact = 0;
caseNameNoBact = [baseName '_NO_BACT'];
problemNoBact = packSimulationProblem(state0NoBact, modelNoBact, schedule, caseNameNoBact, 'NonLinearSolver', nls);
simulatePackedProblem(problemNoBact);

%% Get and compare results
[wsWithClog, statesWithClog] = getPackedSimulatorOutput(problemWithClogging);
[wsNoClog, statesNoClog] = getPackedSimulatorOutput(problemNoClogging);
[wsNoBact, statesNoBact] = getPackedSimulatorOutput(problemNoBact);

% Plot states
figure; plotToolbar(model.G, statesWithClog); title('With Clogging Effects');
figure; plotToolbar(model.G, statesNoClog); title('Without Clogging Effects');
figure; plotToolbar(model.G, statesNoBact); title('Without Bacterial Effects');

% Plot well outputs
figure; plotWellSols({wsWithClog, wsNoClog, wsNoBact});

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