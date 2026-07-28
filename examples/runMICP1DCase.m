%% Basic workflow example for the MICP model
%
% This example demonstrates the complete workflow for creating, running,
% and analyzing a one-dimensional flow system using the MICP mathematical
% model in MATLAB and GNU Octave.
%
% For details on the MICP model, see:
%
% Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E., 2021.
% Practical approaches to study microbially induced calcite precipitation
% at the field scale. International Journal of Greenhouse Gas Control 106,
% 103256.
%
% https://doi.org/10.1016/j.ijggc.2021.103256

% Required modules
mrstModule add ad-blackoil ad-core ad-micp

%% Reservoir geometry, properties, and model parameters
%
% The domain has a length of `domainLength` meters. It is discretized using
% fine cells within the treatment region, except for the first cell, which
% contains the injection well. Beyond the treatment region, the grid is
% gradually coarsened toward the right boundary. The complete domain has
% dimensions `domainLength` x 1 x 1 meters.

% Grid parameters
domainLength = 75;       % Aquifer length, m
treatmentLength = 25;    % Region where MICP is most relevant, m
wellCellSize = 0.5;      % Size of the injection-well cell, m
fineCellSize = 0.25;     % Cell size within the treatment region, m

% The exponential-grid parameters control the coarsening outside the
% treatment region and may need adjustment if the domain parameters change.
xCoordinates = [ ...
    0, ...
    wellCellSize, ...
    wellCellSize + fineCellSize : fineCellSize : treatmentLength, ...
    domainLength * exp(-1.075 : 0.025 : 0)];

xCoordinates(end) = domainLength;

G = tensorGrid(xCoordinates, [0 1], [0 1]);
G = computeGeometry(G);

cellTemplate = ones(G.cells.num, 1);

% Rock properties
initialPermeability = 1e-12 * cellTemplate;  % Permeability, m^2
initialPorosity = 0.2;                       % Porosity, [-]

rock = makeRock( ...
    G, initialPermeability, initialPorosity);

% Fluid properties
fluid.muw = 2.535e-4;          % Water viscosity, Pa s
fluid.bW = @(pressure) ...
    0 * pressure + 1;          % Water formation volume factor, [-]
fluid.rhoWS = 1045;            % Water density, kg/m^3

% MICP model parameters
fluid.rho_b = 35;              % Biofilm density, kg/m^3
fluid.rho_c = 2710;            % Calcite density, kg/m^3
fluid.k_str = 2.6e-10;         % Detachment coefficient, m/(Pa s)
fluid.diffm = 2.1e-9;          % Microorganism diffusion, m^2/s
fluid.diffo = 2.32e-9;         % Oxygen diffusion, m^2/s
fluid.diffu = 1.38e-9;         % Urea diffusion, m^2/s
fluid.alphaL = 1e-3;           % Longitudinal dispersivity, m
fluid.alphaT = 4e-4;           % Transverse dispersivity, m
fluid.eta = 3;                 % Permeability fitting factor, [-]
fluid.k_o = 2e-5;              % Oxygen half-velocity constant, kg/m^3
fluid.k_u = 21.3;              % Urea half-velocity constant, kg/m^3
fluid.mu = 4.17e-5;            % Maximum specific growth rate, 1/s
fluid.mu_u = 0.0161;           % Maximum urease utilization rate, 1/s
fluid.k_a = 8.51e-7;           % Microorganism attachment rate, 1/s
fluid.k_d = 3.18e-7;           % Microorganism death rate, 1/s
fluid.Y = 0.5;                 % Growth yield coefficient, [-]
fluid.Yuc = 1.67;              % Calcite-to-urea yield coefficient, [-]
fluid.F = 0.5;                 % Oxygen consumption factor, [-]
fluid.crit = 0.1;              % Critical porosity, [-]
fluid.kmin = 1e-20;            % Minimum permeability, m^2

% Porosity-permeability relationship
normalizedPorosity = @(porosity) max( ...
    (porosity - fluid.crit) ./ ...
    (initialPorosity - fluid.crit), 0);

fluid.K = @(porosity) max( ...
    (initialPermeability .* ...
    normalizedPorosity(porosity) .^ fluid.eta + ...
    fluid.kmin) .* ...
    initialPermeability ./ ...
    (initialPermeability + fluid.kmin), ...
    fluid.kmin);

% The current MICP implementation considers single-phase water flow but
% derives from an oil-water model. These compatibility fields are therefore
% required by the parent model.
fluid.bO = fluid.bW;
fluid.rhoOS = fluid.rhoWS;

%% Injection strategy, well, outflow boundary, and schedule
%
% Each row of `injectionStrategy` contains:
%
%   1. Phase duration, s
%   2. Simulation timestep, s
%   3. Injection rate, m^3/s
%   4. Microorganism concentration, kg/m^3
%   5. Oxygen concentration, kg/m^3
%   6. Urea concentration, kg/m^3
%
% Different timesteps are used when the well is active and shut in. The
% strategy was selected to produce permeability reduction approximately
% ten meters from the injection well. Optimization techniques could be
% used to design a more efficient strategy.

activeTimeStep = 20 * minute;
shutInTimeStep = hour;

injectionStrategy = [ ...
    20  * hour, activeTimeStep, 2 / day, 0.01, 0,    0; ...
    20  * hour, activeTimeStep, 2 / day, 0,    0,    0; ...
    100 * hour, shutInTimeStep, 0,       0,    0,    0; ...
    20  * hour, activeTimeStep, 2 / day, 0,    0.04, 0; ...
    20  * hour, activeTimeStep, 2 / day, 0,    0,    0; ...
    30  * hour, shutInTimeStep, 0,       0,    0,    0; ...
    20  * hour, activeTimeStep, 2 / day, 0,    0,   60; ...
    20  * hour, activeTimeStep, 2 / day, 0,    0,    0; ...
    30  * hour, shutInTimeStep, 0,       0,    0,    0];

if exist('AD_MICP_TEST', 'var')
    % Use a shorter schedule for regression testing.
    injectionStrategy = [ ...
        20 * hour, activeTimeStep, 2 / day, 0.01, 0.04, 60];
end

numberOfPhases = size(injectionStrategy, 1);

% Create the injection well. Components enter through the leftmost cell.
wellRadius = 0.15;

W = addWell( ...
    [], ...
    G, ...
    rock, ...
    1, ...
    'Type', 'rate', ...
    'Comp_i', [1, 0], ...
    'Val', injectionStrategy(1, 3), ...
    'Radius', wellRadius);

% Add the MICP component controls to the well structure.
W.m = injectionStrategy(1, 4);
W.o = injectionStrategy(1, 5);
W.u = injectionStrategy(1, 6);

% The well is located at the domain boundary. Store the well-cell index so
% the cell-centered velocity can be corrected when computing dispersion
% and detachment.
G.injectionwellonboundary = 1;
G.cellsinjectionwell = 1;

% Apply a constant-pressure outflow boundary on the right side.
boundaryFaceIndices = boundaryFaces(G);

outflowFaces = boundaryFaceIndices( ...
    abs(G.faces.normals(boundaryFaceIndices, 1)) > eps & ...
    G.faces.centroids(boundaryFaceIndices, 1) > ...
    xCoordinates(end - 1));

bc = addBC( ...
    [], outflowFaces, ...
    'pressure', atm, ...
    'sat', [1, 0]);

numberOfBoundaryFaces = size(bc.sat, 1);

bc.m = zeros(numberOfBoundaryFaces, 1);
bc.o = zeros(numberOfBoundaryFaces, 1);
bc.u = zeros(numberOfBoundaryFaces, 1);
bc.b = zeros(numberOfBoundaryFaces, 1);
bc.c = zeros(numberOfBoundaryFaces, 1);

%% Construct simulation schedule

stepsPerPhase = ...
    injectionStrategy(:, 1) ./ injectionStrategy(:, 2);

stepsPerPhase = round(stepsPerPhase);
totalNumberOfSteps = sum(stepsPerPhase);

timesteps = zeros(totalNumberOfSteps, 1);
controlIndices = zeros(totalNumberOfSteps, 1);

firstStep = 1;

for phaseIndex = 1 : numberOfPhases
    lastStep = firstStep + stepsPerPhase(phaseIndex) - 1;

    timesteps(firstStep : lastStep) = ...
        injectionStrategy(phaseIndex, 2);

    controlIndices(firstStep : lastStep) = phaseIndex;

    firstStep = lastStep + 1;
end

schedule = simpleSchedule( ...
    timesteps, 'W', W, 'bc', bc);

baseControl = schedule.control(1);
schedule.control = repmat( ...
    baseControl, numberOfPhases, 1);

for phaseIndex = 1 : numberOfPhases
    schedule.control(phaseIndex).W.val = ...
        injectionStrategy(phaseIndex, 3);

    schedule.control(phaseIndex).W.m = ...
        injectionStrategy(phaseIndex, 4);

    schedule.control(phaseIndex).W.o = ...
        injectionStrategy(phaseIndex, 5);

    schedule.control(phaseIndex).W.u = ...
        injectionStrategy(phaseIndex, 6);
end

schedule.step.control = controlIndices;
schedule.step.val = timesteps;

% Store maximum injected concentrations for state limiting and dynamic
% visualization.
fluid.Comax = max(injectionStrategy(:, 5));
fluid.Cumax = max(injectionStrategy(:, 6));

%% Set up simulation model and solver
%
% The MICP implementation derives from MRST oil-water model. Model and
% solver tolerances can be adjusted after constructing `MICPModel`.

model = MICPModel(G, rock, fluid);
model.toleranceMB = 1e-14;
model.nonlinearTolerance = 1e-14;

solver = getNonLinearSolver(model);
solver.LinearSolver.tolerance = 1e-14;

%% Define initial state
%
% Initially, the domain contains water without dissolved or solid MICP
% components. Different initial conditions, such as an initial biofilm
% fraction, can be assigned here.

state0 = initState(G, W, atm, [1, 0]);

state0.m = zeros(G.cells.num, 1);
state0.o = zeros(G.cells.num, 1);
state0.u = zeros(G.cells.num, 1);
state0.b = zeros(G.cells.num, 1);
state0.c = zeros(G.cells.num, 1);

%% Run the simulation
%
% In MATLAB, `getPlotAfterStepMICP` dynamically visualizes the solution.
% Dynamic plotting is disabled in GNU Octave.

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    [~, states] = simulateScheduleAD( ...
        state0, ...
        model, ...
        schedule, ...
        'NonLinearSolver', solver);
else
    % The last two arguments define the azimuth and elevation angles.
    afterStepFunction = ...
        getPlotAfterStepMICP(state0, model, 0, 270);

    [~, states] = simulateScheduleAD( ...
        state0, ...
        model, ...
        schedule, ...
        'NonLinearSolver', solver, ...
        'afterStepFn', afterStepFunction);
end

%% Process and visualize the results
%
% Results are exported to VTK files for visualization in ParaView. MATLAB
% additionally uses `plotToolbar` for interactive visualization.

outputDirectory = fullfile(pwd, 'vtk_1DCase');

if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end

outputName = fullfile(outputDirectory, 'states');

% Export the results without changing the current working directory.
mrsttovtk(G, states, outputName, '%f');

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    fprintf( ...
        'VTK results written to:\n  %s.pvd\n', ...
        outputName);
    return
end

% Interactive visualization in MATLAB
mrstModule add mrst-gui

figure;
plotToolbar( ...
    G, ...
    states, ...
    'field', 's:1', ...
    'lockCaxis', true);

view([-10, 14]);
axis tight;
colorbar;
caxis([0 1]);

fprintf( ...
    'VTK results written to:\n  %s.pvd\n', ...
    outputName);

%% Copyright notice
%{
Copyright 2021-2026, NORCE Research AS, Computational Geosciences and
Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file. If not, see <http://www.gnu.org/licenses/>.
%}
