%% MRST Simulation for Hydrogen Storage with Bacterial Growth Model
% Description: This script models underground hydrogen storage with microbial 
% activity in a 3D porous medium using compositional fluid properties.
%This case is based on the article:Khoshnevis, N. and Hogeweg, S. and Goncalves Machado, C. and Hagemann, B..
% Numerical Modeling of Bio-Reactive Transport During Underground Hydrogen Storage – a Benchmark Study},
%European Association of Geoscientists and Engineers,vol1,1--5,issn={2214-4609}, 2023.
%https://doi.org/10.3997/2214-4609.202321087
%
% Features:
% - Two-phase system (liquid 'O' and gas 'G')
% - Four components (H2O, H2, CO2, CH4)
% - Microbial activity of archaea
% - Cyclic injection/production schedule
% - Based on EAGE 2023 Benchmark case
%
% Author: [Stéphanie Delage Santacreu]
% Date: [16/09/2025]
% Organization: [Université de Pau et des Pays de l'Adour, E2S UPPA, CNRS, LFCR, UMR5150, Pau, France]

%% Initialization
clear; clc;
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui 

% Simulation parameters
gravity reset on 
nobact = false;         % Set true to disable bacterial effects

%% ============ Grid and Rock Properties =====================
% Define grid dimensions
[nx, ny, nz] = deal(21, 21, 8);       % Grid cells in x, y, z directions
[Lx, Ly, Lz] = deal(1525, 1525, 50);  % Physical dimensions (meters)

% Create grid and shift vertically
G = cartGrid([nx, ny, nz], [Lx, Ly, Lz]);
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + 1000; % Reservoir depth = 1000m
G = computeGeometry(G);

% Define rock properties
perm = [100, 100, 10] .* milli*darcy;  % Permeability in x,y,z
rock = makeRock(G, perm, 0.2);          % Porosity = 0.2

%% Fluid Properties Initialization
% Define compositional fluid
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, ...
                                      {'H2O', 'H2', 'CO2', 'C1'});

% Fluid properties
[rhow, rhog] = deal(999.7 * kilogram/meter^3, 1.2243 * kilogram/meter^3); % Densities
[viscow, viscog] = deal(1.3059 * centi*poise, 0.01763 * centi*poise);     % Viscosities
[cfw, cfg] = deal(5.0015e-5/barsa, 1.0009/barsa);                         % Compressibilities

% Relative permeability parameters
[srw, src] = deal(0.2, 0.05);  % Residual saturations
P0 = 100 * barsa;              % Reference pressure

% Initialize fluid model
fluid = initSimpleADIFluid('phases', 'OG', ...
                           'mu', [viscow, viscog], ...
                           'rho', [rhow, rhog], ...
                           'pRef', P0, ...
                           'c', [cfw, cfg], ...
                           'n', [2, 2], ...
                           'smin', [srw, src]);

% Capillary pressure function
Pe = 0.1 * barsa;
pcOG = @(sw) Pe * sw.^(-1/2);
fluid.pcOG = @(sg) pcOG(max((1 - sg - srw) / (1 - srw), 1e-5));

%% Well Configuration
% Initialize well structures
W1 = []; W2 = []; W3 = []; W5 = [];

% Well placement
n1 = floor(0.5*nx) + 1; 
n2 = floor(0.5*nx) + 1;

% Create wells with different purposes
% 1. CO2-rich injection well
W1 = verticalWell(W1, G, rock, n1, n2, 1:nz-1, ...
                 'comp_i', [0, 1], 'Radius', 0.5, ...
                 'name', 'Injector', 'type', 'rate', ...
                 'Val', 1e6*meter^3/day, 'sign', 1);
W1(1).components = [0.0, 0.6, 0.4, 0.0];  % [H2O, H2, CO2, C1]

% 2. Rest period well (shut-in)
W2 = verticalWell(W2, G, rock, n1, n2, 1:nz-1, ...
                 'compi', [0, 1], 'Radius', 0.5, ...
                 'name', 'Rest', 'type', 'rate', ...
                 'Val', 0.0, 'sign', 1);
W2(1).components = [0.0, 0.95, 0.05, 0.0];

% 3. H2-rich injection well
W3 = verticalWell(W3, G, rock, n1, n2, 1:nz-1, ...
                 'comp_i', [0, 1], 'Radius', 0.5, ...
                 'name', 'Injector', 'type', 'rate', ...
                 'Val', 1e6*meter^3/day, 'sign', 1);
W3(1).components = [0.0, 0.95, 0.05, 0.0];

% 4. Production well
W5 = verticalWell(W5, G, rock, n1, n2, 1:nz-1, ...
                 'compi', [0, 1], 'Radius', 0.5, ...
                 'name', 'Prod', 'type', 'rate', ...
                 'Val', -1e6*meter^3/day, 'sign', -1);
W5(1).components = [0.0, 0.95, 0.05, 0.0];

%% Time Stepping and Schedule
% Define cyclic schedule parameters
ncycles = 6; 
deltaT = 1*day;
nbj_buildUp = 60*day;
nbj_rest = 20*day;
nbj_inject = 30*day;
nbj_idle = 20*day;
nbj_prod = 30*day;
nbj_idle1 = 20*day;

% Create schedule
[schedule, TotalTime, nbuildUp, nrest, ninject, nidle, nprod, nidle1] = ...
    createCyclicScenario2(deltaT, ncycles, nbj_buildUp, nbj_rest, ...
                         nbj_inject, nbj_idle, nbj_prod, nbj_idle1, [W1; W2; W3; W5]);

%% Model Setup
% Equation of state
eosname = 'sw';
model.EOSModel = SoreideWhitsonEquationOfStateModel(G, compFluid, eosname);

% Select backend
if nobact
    backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true, 'rowMajor', true);
else
    backend = DiagonalAutoDiffBackend('modifyOperators', true);
end

% Configure model arguments
arg = {G, rock, fluid, compFluid, true, backend, ...
       'water', false, 'oil', true, 'gas', true, ...
       'bacteriamodel', ~nobact,...
       'liquidPhase', 'O', 'vaporPhase', 'G'};

% Create model
model = BiochemistryModel(arg{:});
model.outputFluxes = false;
model.EOSModel.msalt = 0;

% Setup solvers
lsolve = selectLinearSolverAD(model);
nls = NonLinearSolver();
nls.LinearSolver = lsolve;

%% Initial Conditions
T0 = 40 + 273.15;               % Initial temperature (K)
s0 = [0.2, 0.8];                % Initial saturations (Sw, Sg)
z0 = [0.7, 0.0, 0.02, 0.28];    % Initial composition [H2O, H2, CO2, C1]

% Initialize state
if model.bacteriamodel
    nbact0 = 1; 
    model.nbactMax = 1e8;
    state0 = initCompositionalStateBacteria(model, P0, T0, s0, z0, nbact0, model.EOSModel);
else
    state0 = initCompositionalState(model, P0, T0, s0, z0);
end

%% Run Simulation
% Pack and run simulation problem
name_nbs0 = 'Benchmark2023AEGE_180_pack_NOBACT_6cycles_msalt0';
problem_nbs0 = packSimulationProblem(state0, model, schedule, name_nbs0, 'NonLinearSolver', nls);

if nobact
    % Run without bacteria
    simulatePackedProblem(problem_nbs0, 'restartStep',1);
    [ws_nbs0, states_nbs0] = getPackedSimulatorOutput(problem_nbs0);
else
    % Run with bacteria
    name = 'Benchmark2023AEGE_180_pack_1cycle_n0_1e8_msalt0';
    problem_bs0 = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls); 
    simulatePackedProblem(problem_bs0, 'restartStep',1);
    simulatePackedProblem(problem_nbs0);
    
    % Get results
    [ws_bs0, states_bs0] = getPackedSimulatorOutput(problem_bs0);
    [ws_nbs0, states_nbs0] = getPackedSimulatorOutput(problem_nbs0);
end

%% Post-processing
% Get component indices
namecp = model.EOSModel.getComponentNames();
ncomp = model.EOSModel.getNumberOfComponents;
nT = numel(states_nbs0);
indH2 = find(strcmp(namecp, 'H2'));
indCO2 = find(strcmp(namecp, 'CO2'));
indC1 = find(strcmp(namecp, 'C1'));

% Initialize result arrays
resultFields = {'xH2', 'yH2', 'xCO2', 'yCO2', 'yC1', 'pressure', ...
                'H2_well', 'CO2_well', 'C1_well', 'totMassH2', ...
                'totMassCO2', 'totMassC1', 'FractionMassH2', ...
                'FractionMassCO2', 'FractionMassC1', 'totMassComp'};

for field = resultFields
    eval([field{1} '_nbs0 = zeros(nT, 1);']);
    if ~nobact
        eval([field{1} '_bs0 = zeros(nT, 1);']);
    end
end

% Process results
for i = 1:nT
    xH2_nbs0(i) = max(states_nbs0{i}.x(:, indH2));
    yH2_nbs0(i) = max(states_nbs0{i}.y(:, indH2));
    yC1_nbs0(i) = max(states_nbs0{i}.y(:, indC1));
    xCO2_nbs0(i) = max(states_nbs0{i}.x(:, indCO2));
    yCO2_nbs0(i) = max(states_nbs0{i}.y(:, indCO2));
    pressure_nbs0(i) = mean(states_nbs0{i}.pressure(:));
    H2_well_nbs0(i) = ws_nbs0{i}.H2;
    CO2_well_nbs0(i) = ws_nbs0{i}.CO2;
    C1_well_nbs0(i) = ws_nbs0{i}.C1;
    
    if ~nobact
        xH2_bs0(i) = max(states_bs0{i}.x(:, indH2));
        yH2_bs0(i) = max(states_bs0{i}.y(:, indH2));
        yC1_bs0(i) = max(states_bs0{i}.y(:, indC1));
        xCO2_bs0(i) = max(states_bs0{i}.x(:, indCO2));
        yCO2_bs0(i) = max(states_bs0{i}.y(:, indCO2));
        pressure_bs0(i) = mean(states_bs0{i}.pressure(:));
        H2_well_bs0(i) = ws_bs0{i}.H2;
        CO2_well_bs0(i) = ws_bs0{i}.CO2;
        C1_well_bs0(i) = ws_bs0{i}.C1;
    end
end

% Calculate mass balances
for i = 1:nT
    for j = 1:ncomp
        totMassComp_nbs0(i) = totMassComp_nbs0(i) + ...
                              sum(states_nbs0{i}.FlowProps.ComponentTotalMass{j});
    end
    
    totMassH2_nbs0(i) = sum(states_nbs0{i}.FlowProps.ComponentTotalMass{indH2});
    totMassCO2_nbs0(i) = sum(states_nbs0{i}.FlowProps.ComponentTotalMass{indCO2});
    totMassC1_nbs0(i) = sum(states_nbs0{i}.FlowProps.ComponentTotalMass{indC1});
    
    FractionMassH2_nbs0(i) = totMassH2_nbs0(i) / totMassComp_nbs0(i);
    FractionMassCO2_nbs0(i) = totMassCO2_nbs0(i) / totMassComp_nbs0(i);
    FractionMassC1_nbs0(i) = totMassC1_nbs0(i) / totMassComp_nbs0(i);
    
    if ~nobact
        for j = 1:ncomp
            totMassComp_bs0(i) = totMassComp_bs0(i) + ...
                                 sum(states_bs0{i}.FlowProps.ComponentTotalMass{j});
        end
        
        totMassH2_bs0(i) = sum(states_bs0{i}.FlowProps.ComponentTotalMass{indH2});
        totMassCO2_bs0(i) = sum(states_bs0{i}.FlowProps.ComponentTotalMass{indCO2});
        totMassC1_bs0(i) = sum(states_bs0{i}.FlowProps.ComponentTotalMass{indC1});
        
        FractionMassH2_bs0(i) = totMassH2_bs0(i) / totMassComp_bs0(i);
        FractionMassCO2_bs0(i) = totMassCO2_bs0(i) / totMassComp_bs0(i);
        FractionMassC1_bs0(i) = totMassC1_bs0(i) / totMassComp_bs0(i);
    end
end

%% Calculate Efficiency Metrics
% H2 production efficiency (no bacteria case)
ndeb0 = nbuildUp + nrest;
njcycle = ninject + nidle + nprod + nidle1;

mH2_well_injected_nbs0 = sum(H2_well_nbs0(1:nbuildUp));
mH2_well_produced_nbs0 = 0.0;

for cycle = 1:ncycles
    ndebi = ndeb0 + (cycle-1)*njcycle; 
    nj1 = ndebi + ninject;
    
    mH2_well_injected_nbs0 = mH2_well_injected_nbs0 + ...
                             sum(H2_well_nbs0(ndebi+1:nj1));
    mH2_well_produced_nbs0 = mH2_well_produced_nbs0 + ...
                             sum(H2_well_nbs0(nj1+nidle+1:nj1+nidle+nprod));
end

Efficiency_H2_well_nbs0 = abs(mH2_well_produced_nbs0 ./ mH2_well_injected_nbs0) * 100;
fprintf('In the well, H2 Production Efficiency_nobact, msalt=0 : %.2f%%\n', Efficiency_H2_well_nbs0);

if ~nobact
    % H2 production efficiency (with bacteria case)
    mH2_well_injected_bs0 = sum(H2_well_bs0(1:nbuildUp));
    mH2_well_produced_bs0 = 0.0;
    
    for cycle = 1:ncycles
        ndebi = ndeb0 + (cycle-1)*njcycle; 
        nj1 = ndebi + ninject;
        
        mH2_well_injected_bs0 = mH2_well_injected_bs0 + ...
                                sum(H2_well_bs0(ndebi+1:nj1));
        mH2_well_produced_bs0 = mH2_well_produced_bs0 + ...
                                sum(H2_well_bs0(nj1+nidle+1:nj1+nidle+nprod));
    end
    
    Efficiency_H2_well_bs0 = abs(mH2_well_produced_bs0 ./ mH2_well_injected_bs0) * 100;
    fprintf('In the well, H2 Production Efficiency_bact, msalt=0 : %.2f%%\n', Efficiency_H2_well_bs0);
    
    % Calculate component changes due to bacteria
    H2_loss_percentage_nbs0 = (abs(totMassH2_nbs0 - totMassH2_bs0) ./ totMassH2_nbs0) * 100;
    CO2_loss_percentage_nbs0 = (abs(totMassCO2_nbs0 - totMassCO2_bs0) ./ totMassCO2_nbs0) * 100;
    C1_loss_percentage_nbs0 = (abs(totMassC1_nbs0 - totMassC1_bs0) ./ totMassC1_nbs0) * 100;
    
    fprintf('Total H2 loss due to bacterial effects, msalt=0: %.2f%%\n', H2_loss_percentage_nbs0(end));
    fprintf('Total CO2 loss due to bacterial effects, msalt=0: %.2f%%\n', CO2_loss_percentage_nbs0(end));
    fprintf('Total C1 production due to bacterial effects, msalt=0: %.2f%%\n', C1_loss_percentage_nbs0(end));

    %% Plot Results
    plot_benchmarck_2023AEGE(nT,pressure_nbs0,pressure_bs0,H2_loss_percentage_nbs0,...
    totMassH2_nbs0,totMassH2_bs0,G,states_bs0,nbact0)
end



