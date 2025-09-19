%% Illustrative Case: Efficiency of Structural Trap for Hydrogen Storage
% 
% This example simulates the storage of hydrogen in a saline aquifer and 
% investigates the influence of physical properties such as H2 solubility, 
% and salinity on various factors, including the 
% hydrogen plume, reservoir pressure, caprock tightness, and the recoverability 
% and loss of hydrogen. The model features a 2D dome-shaped aquifer with 
% integrated caprock and bedrock layers, extending horizontally for 50 m 
% and vertically for 50 m.  It represents the first case discussed 
% in the paper: "Phase behavior and black-oil simulations of hydrogen 
% storage in saline aquifers" (Ahmed, E., et al., 2024).
%
% In this simulation, various scenarios are explored regarding the 
% behavior of hydrogen in saline aquifers. The reader is encouraged 
% to examine different input files for specific case studies:
% - **H2STORAGE_RS.DATA**: Models only dissolved gas in brine, with 
%   evaporated water interactions disabled.
% - **H2STORAGE_RS_SALT.DAT**: Represents hydrogen interacting with brine 
%   at a molality of 4 mol/kg.
%
% The upper boundary of the aquifer is determined by the function 
% F(x) = σ + r sin(πx), where σ is set to 25 and r to 5, creating a 
% trapping structure for efficient H2 storage. The containment of injected 
% hydrogen gas relies on sufficiently high entry (capillary) pressure and 
% the low permeability of the caprock and bedrock formations.
%
% Hydrostatic pressure boundary conditions are applied to the lateral 
% boundaries, calculated from a reference pressure p_r at depth z = 0 of 
% 40 bar. No-flux conditions are enforced at the top and bottom boundaries. 
% The initial state of the brine occupies the domain in hydrostatic 
% equilibrium matching the boundary conditions. Physical and fluid parameters 
% are derived from Oko-roafor et al. (2023) and summarized in a related table. 
% The simulation also tests the injection of H2 into both pure water and brine 
% with salinity of 5 mol kg−1.
%--------------------------------------------------------------------------

clearvars; 
mrstModule add biochemistry h2store ad-core ad-blackoil ad-props deckformat mrst-gui upr test-suite spe10

%% Define the case name and read the Eclipse deck file
baseName = 'H2_STORAGE_DOME_TRAP';

%% Use H2STORAGE_RS_SALT.DATA for brine
dataPath = getDatasetPath('h2storage');
dataFile = fullfile(dataPath, 'H2STORAGE_RS.DATA');
deck = readEclipseDeck(dataFile);

%% Notice on Computational Cost
warning('ComputationalCost:High', ...
       ['Please be advised that this is a 10 cycle injection example which often takes some time ', ...
        'to run: reduce cycles for example']);

%% Set up the simulation parameters and model components
[~, ~, state0Bo, modelBo, scheduleBo, ~] = modelForSimple2DAquifer(deck);
bacteriamodel = true;
CLOG = true;

% Convert black-oil model to compositional model
model = convertBlackOilModelToCompositionalModel(modelBo);
state0=convertBlackOilStateToCompositional(modelBo,state0Bo);
state0.components = ensureMinimumFraction(state0.components,6.0e-4);

%% Set up the compositional model with Water and Hydrogen
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, {'H2O', 'H2', 'CO2', 'C1'});
EOS = EquationOfStateModel([], compFluid, 'sw');
T0 = 40 * Kelvin; % Initial temperature in Kelvin
T0 = 273.15 * Kelvin + T0; % Convert to absolute temperature

if (compFluid.getNumberOfComponents>2)
    state0.components(:,3) = state0.components(:,2)/3;
    state0.components(:,4) = state0.components(:,3);
    state0.components(:,2) = state0.components(:,3);
    state0.components(:,1) = 1- sum(state0.components(:,2:end),2);
    state0.components =state0.components.*0+ [0.8480,  1.0e-5,  1.0e-5,   0.1760-0.0240];
    nbact0 = 5; % Initial number of bacteria
    EOS = SoreideWhitsonEquationOfStateModel([], compFluid, 'sw');
    T0 = 44.35 * Kelvin; % Initial temperature in Kelvin
    T0 = 273.15 * Kelvin + T0; % Convert to absolute temperature
end

model.EOSModel = EOS;
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);

if CLOG
    %% Case with Bio-Clogging Effects
    %--------------------------------------------------------------------------
    % Define porosity and permeability reduction due to bacterial growth.
    %--------------------------------------------------------------------------
    % Initial rock properties
    poro0 = model.rock.poro;  % Initial porosity
    perm0 = model.rock.perm(:, 1);  % Initial permeability
    
    % Bacterial concentration parameters
    nc = 120;  % Critical bacterial concentration
    cp = 1.5;  % Clogging coefficient
    scale = 1 + (nbact0 / nc).^2;
    
    % Define porosity multiplier as a function of bacterial concentration
    pvMult_nbact = @(nbact) 1 ./ (1 + (nbact / nc).^2);
    model.fluid.pvMultR = @(p, nbact) pvMult_nbact(nbact);
    poro = @(p, nbact) scale.*poro0 .* pvMult_nbact(nbact);
    model.rock.poro = poro;
    
    % Define permeability as a function of porosity
    tau = @(p, nbact) ((1 - poro0) ./ (1 - poro(p, nbact))).^2 .* (poro(p, nbact) ./ poro0).^3;
    perm = @(p, nbact) perm0 .* tau(p, nbact);
    model.rock.perm = perm;
    
    % Set case name for clogging scenario
    caseNameWithClogging = [baseName '_WITH_CLOGGING'];
end

% Set up additional model properties for biochemistry and flow
arg = {model.G, model.rock, model.fluid, compFluid,...
       false, diagonal_backend, 'oil', true, 'gas', true, ... % Define phases for water-oil system
    'bacteriamodel', bacteriamodel, 'liquidPhase', 'O', ...
    'vaporPhase', 'G'}; % Set phases and EOS model
model = BiochemistryModel(arg{:});
model.gravity = modelBo.gravity;

%% Manually set up well compositions, temperature, and bacteria
schedule = scheduleBo; % Load schedule from the black-oil model
state0.T = state0.T.*0+T0;

% Initialize compositional state with bacteria model
if bacteriamodel
    state0 = initCompositionalStateBacteria(model, state0.pressure, T0, state0.s, state0.components, nbact0, EOS);
else
     state0 = initCompositionalState(model, state0.pressure, T0, state0.s, state0.components, EOS);
end

bc = schedule.control(1).bc; % Get boundary condition from the schedule
cells_bc = sum(model.G.faces.neighbors(bc.face, :), 2);

%% Update schedule controls and boundary conditions
if (compFluid.getNumberOfComponents>2)
    for i = 1:length(schedule.control)
        % We inject 50% H2 and 50% CO2 in the build-up phase
        schedule.control(i).W.compi = [0, 1]; % Well components
        if (strcmp(schedule.control(i).W.name, 'cushion')&&i<11)
            schedule.control(i).W.components = [0.0, 0.1, 0.9 0.0];
        else
            schedule.control(i).W.components = [0.0, 0.95, 0.05, 0.0];
        end
        schedule.control(i).W.T = T0; % Set temperature (not used, but necessary)
        % Update boundary conditions for each control
        schedule.control(i).bc.components = repmat(state0.components(1,:), numel(cells_bc), 1);
        schedule.control(i).bc.sat = repmat(state0.s(1,:), numel(cells_bc), 1);
        schedule.control(i).bc = [];
    end
else
    for i = 1:length(schedule.control)
        % Set component mix for each control, adjusting for specific conditions
        schedule.control(i).W.compi = [0, 1]; % Well components
        if (strcmp(schedule.control(i).W.name, 'cushion'))
            schedule.control(i).W.components = [0.0, 1.0];
        else
            schedule.control(i).W.components = [0.0, 1.0];
        end
        schedule.control(i).W.T = T0; % Set temperature (not used, but necessary)
        % Update boundary conditions for each control
        schedule.control(i).bc.components = repmat(state0.components(1,:), numel(cells_bc), 1); % Set boundary component concentrations
        schedule.control(i).bc.sat = repmat(state0.s(1,:), numel(cells_bc), 1); % Set boundary saturation
    end
end

model.outputFluxes = false;

%% Plot Grid with Wells, Permeability, and Porosity
figure;
subplot(1, 2, 1);
plotCellData(model.G, poro0);
hold on;
plotGrid(model.G, schedule.control(1).W.cells, 'LineStyle', 'none', 'FaceColor', 'red');
title('Porosity');
axis off tight;
subplot(1, 2, 2);
plotCellData(model.G, log10(perm0(:, 1)));
hold on;
plotGrid(model.G, schedule.control(1).W.cells, 'LineStyle', 'none', 'FaceColor', 'red');
title('Permeability (log10)');
axis off tight;
sgtitle('2D Dome-Shaped Aquifer');

%% Initialize the nonlinear solver and select the linear solver
nls = NonLinearSolver();
lsolve = selectLinearSolverAD(model);
nls.LinearSolver = lsolve;

%% Pack and run simulation for case WITH clogging
problemWithClogging = packSimulationProblem(state0, model, schedule, caseNameWithClogging, 'NonLinearSolver', nls);
simulatePackedProblem(problemWithClogging);

%% Set up and run simulation for case WITHOUT clogging
modelNoClogging = model;
modelNoClogging.rock.perm = perm0;  % Restore original permeability
modelNoClogging.rock.poro = poro0;  % Restore original porosity
modelNoClogging.fluid.pvMultR = @(p, nbact) 1;  % Remove clogging effect
bacteriamodel = true;

arg = {modelNoClogging.G, modelNoClogging.rock, modelNoClogging.fluid, compFluid,...
       false, diagonal_backend, 'oil', true, 'gas', true, ...
    'bacteriamodel', bacteriamodel,'liquidPhase', 'O', ...
    'vaporPhase', 'G'};
modelNoClogging = BiochemistryModel(arg{:});

state0NoClogging = state0;
caseNameNoClogging = [baseName '_NO_CLOGGING'];

problemNoClogging = packSimulationProblem(state0NoClogging, modelNoClogging, schedule, caseNameNoClogging, 'NonLinearSolver', nls);
simulatePackedProblem(problemNoClogging, 'restartStep', 1);

%% Set up and run simulation wthout bacteria
modelNoBact = model;
modelNoBact.rock.perm = perm0;  % Restore original permeability
modelNoBact.rock.poro = poro0;  % Restore original porosity
modelNoBact.fluid.pvMultR = @(p, nbact) 1;  % Remove clogging effect
bacteriamodel = false;

arg = {modelNoBact.G, modelNoBact.rock, modelNoBact.fluid, compFluid,...
       false, diagonal_backend, 'oil', true, 'gas', true, ...
    'bacteriamodel', bacteriamodel, 'liquidPhase', 'O', ...
    'vaporPhase', 'G'};
modelNoBact = BiochemistryModel(arg{:});

state0NoBact = state0;
state0NoBact.nbact = 0;  % No bacteria in this case
caseNameNoBact = [baseName '_NO_BACT'];

problemNoBact = packSimulationProblem(state0NoBact, modelNoBact, schedule, caseNameNoBact, 'NonLinearSolver', nls);
simulatePackedProblem(problemNoBact, 'restartStep', 1);
%% Get and compare results from both cases
[wsWithClog, statesWithClog] = getPackedSimulatorOutput(problemWithClogging);
[wsNoClog, statesNoClog] = getPackedSimulatorOutput(problemNoClogging);
[wsNoBact, statesNoBact] = getPackedSimulatorOutput(problemNoBact);


% Plot states comparison
figure;
plotToolbar(model.G, statesWithClog);
title('With Clogging Effects');
figure;
plotToolbar(model.G, statesNoClog);
title('Without Clogging Effects');
figure;
plotToolbar(model.G, statesNoBact);
title('Without Bacterial Effects');

% Plot well output comparison
figure;
plotWellSols({wsWithClog, wsNoClog, wsNoBact});