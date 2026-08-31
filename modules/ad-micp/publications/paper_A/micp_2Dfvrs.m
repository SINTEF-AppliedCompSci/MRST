%% Two-dimensional vertical rectangular MICP treatment
%
% Set up and solve the two-dimensional vertical rectangular-flow system
% (2Dfvrs) using two alternative injection strategies.
%
% In MATLAB, this example produces Figure 8 from:
%
% Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E., 2021. Practical
% approaches to study microbially induced calcite precipitation at the
% field scale. International Journal of Greenhouse Gas Control 106, 103256.
%
% https://doi.org/10.1016/j.ijggc.2021.103256
%
% In GNU Octave, the example writes the results to the
% `vtk_micp_2Dfvrs` directory for visualization in ParaView.
%
% The example assumes that MRST and the `ad-micp` module are available on
% the MATLAB or GNU Octave path.

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

%% Required modules

mrstModule add ad-blackoil ad-core ad-micp

%% Grid
aquiferLength = 500; % Aquifer length, m
aquiferHeight = 30;  % Aquifer height, m
if exist('AD_MICP_TEST', 'var')
    % Smaller case for the test.
    xCoordinates = [ ...
        0 : 10 : 200, ...
        aquiferLength .* exp(-0.9 : 0.1 : 0)];
    xCoordinates = unique(xCoordinates, 'stable');
    xCoordinates(end) = aquiferLength;
    zCoordinates = 0 : 3 : aquiferHeight;
else
    xCoordinates = [ ...
        0 : 0.5 : 200, ...
        aquiferLength .* exp(-0.9 : 0.05 : 0)];
    xCoordinates = unique(xCoordinates, 'stable');
    xCoordinates(end) = aquiferLength;
    zCoordinates = 0 : 0.5 : aquiferHeight;
end

G = tensorGrid(xCoordinates, zCoordinates, [0 1]);
G = computeGeometry(G);

cellCentroids = G.cells.centroids;
cellTemplate = ones(G.cells.num, 1);

%% Rock properties

initialPermeability = 2e-14 .* cellTemplate; % Permeability, m^2
initialPorosity = 0.15;                      % Porosity, [-]

rock = makeRock(G, initialPermeability, initialPorosity);

%% Fluid and MICP model properties

fluid.muw = 2.535e-4;       % Water viscosity, Pa s
fluid.bW = @(pressure) ...
    0 .* pressure + 1;      % Water formation volume factor, [-]
fluid.bO = fluid.bW;        % Compatibility field, [-]
fluid.rhoWS = 1045;         % Water density, kg/m^3
fluid.rhoOS = 479;          % CO2 density, kg/m^3

fluid.rho_b = 35;           % Biofilm density, kg/m^3
fluid.rho_c = 2710;         % Calcite density, kg/m^3
fluid.k_str = 2.6e-10;      % Detachment coefficient, m/(Pa s)
fluid.diffm = 2.1e-9;       % Microorganism diffusion, m^2/s
fluid.diffo = 2.32e-9;      % Oxygen diffusion, m^2/s
fluid.diffu = 1.38e-9;      % Urea diffusion, m^2/s
fluid.alphaL = 1e-3;        % Longitudinal dispersivity, m
fluid.alphaT = 4e-4;        % Transverse dispersivity, m
fluid.eta = 3;              % Permeability fitting factor, [-]
fluid.k_o = 2e-5;           % Oxygen half-velocity constant, kg/m^3
fluid.k_u = 21.3;           % Urea half-velocity constant, kg/m^3
fluid.mu = 4.17e-5;         % Maximum specific growth rate, 1/s
fluid.mu_u = 0.0161;        % Maximum urease utilization rate, 1/s
fluid.k_a = 8.51e-7;        % Microorganism attachment rate, 1/s
fluid.k_d = 3.18e-7;        % Microorganism death rate, 1/s
fluid.Y = 0.5;              % Growth yield coefficient, [-]
fluid.Yuc = 1.67;           % Calcite-to-urea yield coefficient, [-]
fluid.F = 0.5;              % Oxygen consumption factor, [-]
fluid.crit = 0.1;           % Critical porosity, [-]
fluid.kmin = 1e-20;         % Minimum permeability, m^2

% Porosity-permeability relationship
normalizedPorosity = @(currentPorosity) max( ...
    (currentPorosity - fluid.crit) ./ ...
    (initialPorosity - fluid.crit), 0);

fluid.K = @(currentPorosity) max( ...
    (initialPermeability .* ...
    normalizedPorosity(currentPorosity) .^ fluid.eta + ...
    fluid.kmin) .* initialPermeability ./ ...
    (initialPermeability + fluid.kmin), fluid.kmin);

%% Gravity

gravity reset
gravity on
gravity y

gravityMagnitude = norm(gravity);

%% Open boundary condition
%
% Apply a hydrostatic pressure condition at the right boundary.

boundaryFaceIndices = boundaryFaces(G);

outflowFaces = boundaryFaceIndices( ...
    abs(G.faces.normals(boundaryFaceIndices, 1)) > eps & ...
    G.faces.centroids(boundaryFaceIndices, 1) > ...
    aquiferLength - 0.01);

boundaryPressure = ...
    G.faces.centroids(outflowFaces, 2) .* ...
    fluid.rhoWS .* gravityMagnitude;

bc = addBC([], outflowFaces, 'pressure', boundaryPressure, ...
                                      'sat', [0 0]);

numberOfBoundaryFaces = size(bc.sat, 1);

bc.m = zeros(numberOfBoundaryFaces, 1);
bc.o = zeros(numberOfBoundaryFaces, 1);
bc.u = zeros(numberOfBoundaryFaces, 1);
bc.b = zeros(numberOfBoundaryFaces, 1);
bc.c = zeros(numberOfBoundaryFaces, 1);

%% Injection strategy
%
% Each row of `injectionStrategy` contains:
%
%   1. Phase duration, s
%   2. Simulation timestep, s
%   3. Injection rate, m^3/s
%   4. Microorganism concentration, kg/m^3
%   5. Oxygen concentration, kg/m^3
%   6. Urea concentration, kg/m^3

activeTimeStep = hour;

injectionStrategy = [ ...
     15 * hour, activeTimeStep, 5e-3, 0.01, 0,    0; ...
     11 * hour, activeTimeStep, 5e-3, 0,    0,    0; ...
     74 * hour, 2 * hour,       0,    0,    0,    0; ...
     30 * hour, activeTimeStep, 5e-3, 0,    0.04, 0; ...
      5 * hour, activeTimeStep, 5e-3, 0,    0,    0; ...
     25 * hour, 5 * hour,       0,    0,    0,    0; ...
     40 * hour, activeTimeStep, 5e-3, 0,    0,  300; ...
     10 * hour, activeTimeStep, 5e-3, 0,    0,    0; ...
    390 * hour, 26 * hour,      0,    0,    0,    0];

if exist('AD_MICP_TEST', 'var')
    injectionStrategy = [ ...
        2 * hour, activeTimeStep, 5e-3, 0.01, 0.04, 300];
end

numberOfPhases = size(injectionStrategy, 1);

% Maximum injected oxygen and urea concentrations
fluid.Comax = max(injectionStrategy(:, 5));
fluid.Cumax = max(injectionStrategy(:, 6));

%% Construct timestep and control arrays

stepsPerPhase = round( ...
    injectionStrategy(:, 1) ./ injectionStrategy(:, 2));

totalNumberOfSteps = sum(stepsPerPhase);

timesteps = zeros(totalNumberOfSteps, 1);
controlIndices = zeros(totalNumberOfSteps, 1);

firstStep = 1;

for phaseIndex = 1 : numberOfPhases
    lastStep = firstStep + stepsPerPhase(phaseIndex) - 1;

    timesteps(firstStep : lastStep) = ...
        injectionStrategy(phaseIndex, 2);

    controlIndices(firstStep : lastStep) = ...
        phaseIndex;

    firstStep = lastStep + 1;
end

%% Strategy A
%
% Inject all water and MICP components through the lower section of the
% well.

wellRadius = 0.15;
allCellIndices = (1 : G.cells.num)';

strategyACells = allCellIndices( ...
    cellCentroids(:, 1) < 1.1 .* cellCentroids(1, 1) & ...
    cellCentroids(:, 2) < aquiferHeight / 10);

wellA = addWell([], G, rock, strategyACells, 'Type', 'rate', ...
    'Comp_i', [1, 0], 'Val', injectionStrategy(1, 3), ...
    'Radius', wellRadius, 'dir', 'y');

wellA.m = injectionStrategy(1, 4);
wellA.o = injectionStrategy(1, 5);
wellA.u = injectionStrategy(1, 6);

G.injectionwellonboundary = 1;
G.cellsinjectionwell = strategyACells;

scheduleA = simpleSchedule(timesteps, 'W', wellA, 'bc', bc);

baseControlA = scheduleA.control(1);
scheduleA.control = repmat(baseControlA, numberOfPhases, 1);

for phaseIndex = 1 : numberOfPhases
    scheduleA.control(phaseIndex).W.val = ...
        injectionStrategy(phaseIndex, 3);

    scheduleA.control(phaseIndex).W.m = ...
        injectionStrategy(phaseIndex, 4);

    scheduleA.control(phaseIndex).W.o = ...
        injectionStrategy(phaseIndex, 5);

    scheduleA.control(phaseIndex).W.u = ...
        injectionStrategy(phaseIndex, 6);
end

scheduleA.step.control = controlIndices;
scheduleA.step.val = timesteps;

%% Create model and initial state for strategy A

modelA = MICPModel(G, rock, fluid);
modelA.toleranceMB = 1e-14;
modelA.nonlinearTolerance = 1e-12;

initialPressure = ...
    cellCentroids(:, 2) .* fluid.rhoWS .* gravityMagnitude;

state0A = initState(G, wellA, initialPressure, [1, 0]);

state0A.m = zeros(G.cells.num, 1);
state0A.o = zeros(G.cells.num, 1);
state0A.u = zeros(G.cells.num, 1);
state0A.b = zeros(G.cells.num, 1);
state0A.c = zeros(G.cells.num, 1);

%% Simulate strategy A

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    [~, statesA] = simulateScheduleAD( ...
        state0A, modelA, scheduleA);
else
    afterStepFunction = ...
        getPlotAfterStepMICP(state0A, modelA, 0, 270);

    [~, statesA] = simulateScheduleAD( ...
        state0A, modelA, scheduleA, ...
        'afterStepFn', afterStepFunction);
end

finalStateA = statesA{end};

%% Strategy B
%
% Split the well into an upper and lower completion. MICP components are
% injected through the upper completion, while the lower completion
% injects water only.

upperWellFraction = 1 / 10;
lowerWellFraction = 1 - upperWellFraction;

upperStrategyBCells = allCellIndices( ...
    cellCentroids(:, 1) < 1.1 .* cellCentroids(1, 1) & ...
    cellCentroids(:, 2) < upperWellFraction .* aquiferHeight);

lowerStrategyBCells = allCellIndices( ...
    cellCentroids(:, 1) < 1.1 .* cellCentroids(1, 1) & ...
    cellCentroids(:, 2) > upperWellFraction .* aquiferHeight);

wellB = addWell([], G, rock, upperStrategyBCells, 'Type', 'rate', ...
    'Comp_i', [1, 0], 'Val', upperWellFraction .* ...
    injectionStrategy(1, 3), 'Radius', wellRadius, 'dir', 'y');

wellB = addWell(wellB, G, rock, lowerStrategyBCells, 'Type', 'rate', ...
    'Comp_i', [1, 0], 'Val', lowerWellFraction .* ...
    injectionStrategy(1, 3), 'Radius', wellRadius, 'dir', 'y');

for wellIndex = 1 : numel(wellB)
    wellB(wellIndex).m = 0;
    wellB(wellIndex).o = 0;
    wellB(wellIndex).u = 0;
end

wellB(1).m = injectionStrategy(1, 4);
wellB(1).o = injectionStrategy(1, 5);
wellB(1).u = injectionStrategy(1, 6);

G.injectionwellonboundary = 1;
G.cellsinjectionwell = [upperStrategyBCells; lowerStrategyBCells];

scheduleB = simpleSchedule(timesteps, 'W', wellB, 'bc', bc);

baseControlB = scheduleB.control(1);
scheduleB.control = repmat(baseControlB, numberOfPhases, 1);

for phaseIndex = 1 : numberOfPhases
    scheduleB.control(phaseIndex).W(1).val = ...
        upperWellFraction .* injectionStrategy(phaseIndex, 3);

    scheduleB.control(phaseIndex).W(2).val = ...
        lowerWellFraction .* injectionStrategy(phaseIndex, 3);

    scheduleB.control(phaseIndex).W(1).m = ...
        injectionStrategy(phaseIndex, 4);

    scheduleB.control(phaseIndex).W(1).o = ...
        injectionStrategy(phaseIndex, 5);

    scheduleB.control(phaseIndex).W(1).u = ...
        injectionStrategy(phaseIndex, 6);
end

scheduleB.step.control = controlIndices;
scheduleB.step.val = timesteps;

%% Create model and initial state for strategy B

modelB = MICPModel(G, rock, fluid);
modelB.toleranceMB = 1e-14;
modelB.nonlinearTolerance = 1e-12;

state0B = initState(G, wellB, initialPressure, [1, 0]);

state0B.m = zeros(G.cells.num, 1);
state0B.o = zeros(G.cells.num, 1);
state0B.u = zeros(G.cells.num, 1);
state0B.b = zeros(G.cells.num, 1);
state0B.c = zeros(G.cells.num, 1);

%% Simulate strategy B

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    [~, statesB] = simulateScheduleAD( ...
        state0B, modelB, scheduleB);
else
    afterStepFunction = ...
        getPlotAfterStepMICP(state0B, modelB, 0, 270);

    [~, statesB] = simulateScheduleAD( ...
        state0B, modelB, scheduleB, ...
        'afterStepFn', afterStepFunction);
end

finalStateB = statesB{end};

%% Export results in GNU Octave

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    outputDirectory = fullfile(pwd, 'vtk_micp_2Dfvrs');

    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    strategyAOutput = fullfile(outputDirectory, 'stateA');
    strategyBOutput = fullfile(outputDirectory, 'stateB');

    mrsttovtk(G, statesA, strategyAOutput, '%f');
    mrsttovtk(G, statesB, strategyBOutput, '%f');

    fprintf(['VTK results written to:\n' ...
             '  %s.pvd\n' ...
             '  %s.pvd\n'], ...
             strategyAOutput, strategyBOutput);

    return
end

%% Figure 8

lW = 3;
fS = 8;

figure;

colorMap = flipud(jet);
numberOfColors = size(colorMap, 1);
colorMap = colorMap(round(70 * numberOfColors / 256) : end, :);

set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6.83 1.85]);

n1 = subplot(1, 2, 1);
view(0, 270);
colormap(n1, colorMap);
caxis([0 100]);
cb = colorbar;
title(cb, '$\%$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.48 0.25 0.01 0.5], 'YTick', [0 25 50 75 100]);
xlabel({'x [m]'; '(a)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');

porosityA = max(initialPorosity - finalStateA.c - finalStateA.b, ...
                modelA.minimumPorosity);
s = plotCellData(G, 100 .* (1 - fluid.K(porosityA) ./ initialPermeability));
s.EdgeColor = 'none';

title('Permeability reduction (using strategy A)', 'FontSize', fS, ...
                              'FontName', 'Arial', 'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', 0 : 100 : aquiferLength, ...
                             'color', 'none', 'FontName', 'Arial');
zlim([0, aquiferHeight]);
line([90 90], [0 aquiferHeight / 10], [1 1], 'Color', '[0 0 0]', ...
                                        'LineStyle', '-', 'LineWidth', lW);
line([110 110], [0 aquiferHeight / 10], [1 1], 'Color', '[0 0 0]', ...
                                        'LineStyle', '-', 'LineWidth', lW);
line([90 110], [aquiferHeight / 10 aquiferHeight / 10], [1 1], ...
                   'Color', '[0 0 0]', 'LineStyle', '-', 'LineWidth', lW);
line([90 110], [0 0], [1 1], 'Color', '[0 0 0]', ...
                                        'LineStyle', '-', 'LineWidth', lW);

n2 = subplot(1, 2, 2);
view(0, 270);
colormap(n2, colorMap);
caxis([0 100]);
cb = colorbar;
title(cb, '$\%$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.93 0.25 0.01 0.5], 'YTick', [0 25 50 75 100]);
xlabel({'x [m]' ; '(b)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');

porosityB = max(initialPorosity - finalStateB.c - finalStateB.b, ...
                modelB.minimumPorosity);
s = plotCellData(G, 100 .* (1 - fluid.K(porosityB) ./ initialPermeability));
s.EdgeColor = 'none';

title('Permeability reduction (using strategy B)', 'FontSize', fS, ...
                              'FontName', 'Arial', 'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', 0 : 100 : aquiferLength, ...
                             'color', 'none', 'FontName', 'Arial');
zlim([0, aquiferHeight]);
line([90 90], [0 aquiferHeight / 10], [1 1], 'Color', '[0 0 0]', ...
                                        'LineStyle', '-', 'LineWidth', lW);
line([110 110], [0 aquiferHeight / 10], [1 1], 'Color', '[0 0 0]', ...
                                        'LineStyle', '-', 'LineWidth', lW);
line([90 110], [aquiferHeight / 10 aquiferHeight / 10], [1 1], ...
                   'Color', '[0 0 0]', 'LineStyle', '-', 'LineWidth', lW);
line([90 110], [0 0], [1 1], 'Color', '[0 0 0]', ...
                                        'LineStyle', '-', 'LineWidth', lW);

% print -depsc2 Fig8.eps
