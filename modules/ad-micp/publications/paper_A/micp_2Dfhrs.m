%% Two-dimensional horizontal rectangular MICP treatment
%
% Set up and solve the two-dimensional horizontal rectangular-flow system
% (2Dfhrs).
%
% In MATLAB, this example produces Figure 7 from:
%
% Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E., 2021. Practical
% approaches to study microbially induced calcite precipitation at the
% field scale. International Journal of Greenhouse Gas Control 106, 103256.
%
% https://doi.org/10.1016/j.ijggc.2021.103256
%
% In GNU Octave, the example writes the results to the
% `vtk_micp_2Dfhrs` directory for visualization in ParaView.
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

distmeshPath = fullfile(ROOTDIR, 'utils', '3rdparty', 'DistMesh');
% Download DistMesh if it is not already available.
if ~isfolder(distmeshPath)
    distmeshUrl = ...
        'https://github.com/popersson/DistMesh/archive/refs/tags/v1.1.zip';
    distmeshArchive = fullfile(tempdir, 'distmesh.zip');
    extractionDirectory = tempname;
    mkdir(extractionDirectory);
    if isfile(distmeshArchive)
        delete(distmeshArchive);
    end
    urlwrite(distmeshUrl, distmeshArchive);
    unzip(distmeshArchive, extractionDirectory);
    extractedDistmeshPath = fullfile( ...
        extractionDirectory, 'DistMesh-1.1');
    if isfolder(distmeshPath)
        rmdir(distmeshPath, 's');
    end
    movefile(extractedDistmeshPath, distmeshPath);
    if isfile(distmeshArchive)
        delete(distmeshArchive);
    end
    if isfolder(extractionDirectory)
        rmdir(extractionDirectory, 's');
    end
end

mrstPath('reregister', 'distmesh', distmeshPath);
mrstModule add ad-blackoil ad-core ad-micp distmesh

%% Grid

aquiferLength = 75;  % Aquifer half-length, m
aquiferWidth = 10;   % Aquifer half-width, m

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    % Fix the random seed to improve mesh reproducibility.
    rand('seed', 0);
    randn('seed', 0);

    if exist('AD_MICP_TEST', 'var')
        % Smaller case for the test.
        transitionRadius = 25; % Grid-refinement transition radius, m
        minimumGridSize = 3;   % Minimum grid size, m
        maximumGridSize = 10;  % Maximum grid size, m
    else
        transitionRadius = 25; % Grid-refinement transition radius, m
        minimumGridSize = 1;   % Minimum grid size, m
        maximumGridSize = 10;  % Maximum grid size, m
    end

    distanceFunction = @(points) drectangle( ...
        points, -aquiferLength, aquiferLength, ...
        -aquiferWidth, aquiferWidth);

    sizeFunction = @(points) ...
        min(minimumGridSize + ...
        0.3 .* abs(dcircle(points, 0, 0, 0)), minimumGridSize) .* ...
        (abs(dcircle(points, 0, 0, 0)) < transitionRadius) + ...
        min(minimumGridSize + ...
        0.3 .* abs(dcircle(points, 0, 0, transitionRadius)), ...
        maximumGridSize) .* ...
        (abs(dcircle(points, 0, 0, 0)) >= transitionRadius);

    boundingBox = [ ...
        -aquiferLength, -aquiferWidth; ...
         aquiferLength,  aquiferWidth];

    fixedPoints = [ ...
        -aquiferLength, -aquiferWidth; ...
         aquiferLength, -aquiferWidth; ...
        -aquiferLength,  aquiferWidth; ...
         aquiferLength,  aquiferWidth; ...
         0,               0];

    if exist('AD_MICP_TEST', 'var')
        % Disable DistMesh figures during automated tests.
        previousFigureVisibility = ...
            get(0, 'defaultFigureVisible');
        set(0, 'defaultFigureVisible', 'off');
        figureVisibilityCleanup = onCleanup(@() ...
            set(0, 'defaultFigureVisible', ...
            previousFigureVisibility));
        existingFigures = findall(0, 'Type', 'figure');
    end

    [gridPoints, gridTriangles] = distmesh2d( ...
        distanceFunction, sizeFunction, minimumGridSize, ...
        boundingBox, fixedPoints);

    G = makeLayeredGrid( ...
        pebi(triangleGrid(gridPoints, gridTriangles)), 1);

    if exist('AD_MICP_TEST', 'var')
        % Delete only figures created by DistMesh.
        newFigures = setdiff( ...
            findall(0, 'Type', 'figure'), ...
            existingFigures);
        if ~isempty(newFigures)
            delete(newFigures);
        end
        % Restore the original default figure visibility immediately.
        clear figureVisibilityCleanup
    else
        close
    end
else
    [coarseX, coarseY] = meshgrid( ...
        -aquiferLength : 5 : aquiferLength, ...
        -aquiferWidth : 5 : aquiferWidth);

    [leftFineX, leftFineY] = meshgrid( ...
        -30 : 0.5 : -1.5, ...
        -aquiferWidth : 0.5 : aquiferWidth);

    [rightFineX, rightFineY] = meshgrid( ...
        1.5 : 0.5 : 30, ...
        -aquiferWidth : 0.5 : aquiferWidth);

    radialCoordinates = -1 : 0.125 : 0;
    wellPoints = [];

    for radius = 4 .* exp(radialCoordinates)
        [xCoordinates, yCoordinates, ~] = cylinder(radius, 50);
        wellPoints = [wellPoints, ...
                      [xCoordinates(1, :); yCoordinates(1, :)]]; %#ok<AGROW>
    end

    wellPoints = [wellPoints, [0; 0]];

    gridPoints = unique([ ...
        wellPoints'; ...
        coarseX(:), coarseY(:); ...
        leftFineX(:), leftFineY(:); ...
        rightFineX(:), rightFineY(:)], ...
        'rows');

    G = triangleGrid(gridPoints);
    G = computeGeometry(G);
    G = makeLayeredGrid(pebi(G), 1);
end

G = computeGeometry(G);

cellCentroids = G.cells.centroids;
cellTemplate = ones(G.cells.num, 1);

%% Rock properties

initialPermeability = 1e-12 .* cellTemplate; % Permeability, m^2
initialPorosity = 0.2;                       % Porosity, [-]

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
%
% The first nine phases comprise the first MICP treatment, and the final
% nine phases repeat the strategy for the second treatment.

activeTimeStep = hour;
shutInTimeStep = 10 * hour;

injectionStrategy = [ ...
     20 * hour, activeTimeStep, 7.2e-4, 0.01, 0,    0; ...
     20 * hour, activeTimeStep, 7.2e-4, 0,    0,    0; ...
    100 * hour, shutInTimeStep, 0,      0,    0,    0; ...
     20 * hour, activeTimeStep, 7.2e-4, 0,    0.04, 0; ...
     20 * hour, activeTimeStep, 7.2e-4, 0,    0,    0; ...
     50 * hour, shutInTimeStep, 0,      0,    0,    0; ...
     20 * hour, activeTimeStep, 7.2e-4, 0,    0,  300; ...
     20 * hour, activeTimeStep, 7.2e-4, 0,    0,    0; ...
    230 * hour, shutInTimeStep, 0,      0,    0,    0; ...
     20 * hour, activeTimeStep, 7.2e-4, 0.01, 0,    0; ...
     20 * hour, activeTimeStep, 7.2e-4, 0,    0,    0; ...
    100 * hour, shutInTimeStep, 0,      0,    0,    0; ...
     20 * hour, activeTimeStep, 7.2e-4, 0,    0.04, 0; ...
     20 * hour, activeTimeStep, 7.2e-4, 0,    0,    0; ...
     50 * hour, shutInTimeStep, 0,      0,    0,    0; ...
     20 * hour, activeTimeStep, 7.2e-4, 0,    0,  300; ...
     20 * hour, activeTimeStep, 7.2e-4, 0,    0,    0; ...
    230 * hour, shutInTimeStep, 0,      0,    0,    0];

if exist('AD_MICP_TEST', 'var')
    injectionStrategy = [ ...
        5 * hour, activeTimeStep, 7.2e-4, 0.01, 0.04, 300];
end

numberOfPhases = size(injectionStrategy, 1);

%% Injection well

wellRadius = 0.15;

[~, wellCell] = min( ...
    cellCentroids(:, 1) .^ 2 + cellCentroids(:, 2) .^ 2);

W = addWell([], G, rock, wellCell, 'Type', 'rate', ...
    'Comp_i', [1, 0], 'Val', injectionStrategy(1, 3), ...
    'Radius', wellRadius);

W.m = injectionStrategy(1, 4);
W.o = injectionStrategy(1, 5);
W.u = injectionStrategy(1, 6);

% The injection well is not located on the model boundary.
G.injectionwellonboundary = 0;

%% Outflow boundary condition

boundaryFaceIndices = boundaryFaces(G);

outflowFaces = boundaryFaceIndices( ...
    abs(G.faces.normals(boundaryFaceIndices, 1)) > eps & ...
    abs(G.faces.centroids(boundaryFaceIndices, 1)) > ...
    aquiferLength - 0.01);

bc = addBC([], outflowFaces, 'pressure', atm, 'sat', [0 0]);

numberOfBoundaryFaces = size(bc.sat, 1);

bc.m = zeros(numberOfBoundaryFaces, 1);
bc.o = zeros(numberOfBoundaryFaces, 1);
bc.u = zeros(numberOfBoundaryFaces, 1);
bc.b = zeros(numberOfBoundaryFaces, 1);
bc.c = zeros(numberOfBoundaryFaces, 1);

%% Construct simulation schedule

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

    controlIndices(firstStep : lastStep) = phaseIndex;

    firstStep = lastStep + 1;
end

schedule = simpleSchedule(timesteps, 'W', W, 'bc', bc);

baseControl = schedule.control(1);
schedule.control = repmat(baseControl, numberOfPhases, 1);

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

phaseEndSteps = cumsum(stepsPerPhase);

% Maximum injected oxygen and urea concentrations
fluid.Comax = max(injectionStrategy(:, 5));
fluid.Cumax = max(injectionStrategy(:, 6));

%% Create model

model = MICPModel(G, rock, fluid);
model.toleranceMB = 1e-14;
model.nonlinearTolerance = 1e-14;

%% Initial state

state0 = initState(G, W, atm, [1, 0]);

state0.m = zeros(G.cells.num, 1);
state0.o = zeros(G.cells.num, 1);
state0.u = zeros(G.cells.num, 1);
state0.b = zeros(G.cells.num, 1);
state0.c = zeros(G.cells.num, 1);

%% Run the simulation

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    [~, states] = simulateScheduleAD(state0, model, schedule);
else
    afterStepFunction = getPlotAfterStepMICP(state0, model, 0, 90);

    [~, states] = simulateScheduleAD(state0, model, schedule, ...
                             'afterStepFn', afterStepFunction);
end

%% Export results in GNU Octave

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    outputDirectory = fullfile(pwd, 'vtk_micp_2Dfhrs');

    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    outputName = fullfile(outputDirectory, 'states_2Dfhrs');
    mrsttovtk(G, states, outputName, '%f');

    fprintf('VTK results written to:\n  %s.pvd\n', outputName);
    return
end

%% Figure 7

lW = 2;
fS = 8;

figure;

colorMap = flipud(jet);
numberOfColors = size(colorMap, 1);
colorMap = colorMap(round(70 * numberOfColors / 256) : end, :);

set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6.83 1.85]);

% In the full example, phase 9 marks the end of the first treatment. In the
% shortened test case, use the final available phase.
firstTreatmentPhase = min(9, numberOfPhases);
firstTreatmentState = phaseEndSteps(firstTreatmentPhase);

n1 = subplot(1, 2, 1);
colormap(n1, colorMap);
caxis([0 100]);
axis equal tight
colorbar()
cb = colorbar;
title(cb, '$\%$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'location', 'northoutside', 'YTick', [0 25 50 75 100]);
xlabel({'x [m]' ; '(a)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('y [m]', 'FontSize', fS, 'FontName', 'Arial');

firstTreatmentPorosity = max( ...
    initialPorosity - states{firstTreatmentState}.c - ...
    states{firstTreatmentState}.b, model.minimumPorosity);

s = plotCellData(G, 100 .* (1 - ...
    fluid.K(firstTreatmentPorosity) ./ initialPermeability));

s.EdgeColor = 'none';
title('Permeability reduction (after phase I)', 'FontSize', fS, ...
                              'FontName', 'Arial', 'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', -aquiferLength : 25 : aquiferLength, ...
    'YTick', -aquiferWidth : 10 : aquiferWidth, 'color', 'none', ...
    'FontName', 'Arial');
ylim([-aquiferWidth, aquiferWidth]);
rectangle('Position', [10, -aquiferWidth, 5, 2 * aquiferWidth], ...
              'LineWidth', lW, 'LineStyle', '-', 'edgecolor', '[0 0 0]');

n2 = subplot(1, 2, 2);
axis equal tight
colormap(n2, colorMap);
caxis([0 100]);
cb = colorbar;
title(cb, '$\%$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'location', 'northoutside', 'YTick', [0 25 50 75 100]);
xlabel({'x [m]' ; '(b)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('y [m]', 'FontSize', fS, 'FontName', 'Arial');

finalPorosity = max( ...
    initialPorosity - states{end}.c - states{end}.b, ...
    model.minimumPorosity);

s = plotCellData(G, 100 .* (1 - ...
    fluid.K(finalPorosity) ./ initialPermeability));

s.EdgeColor = 'none';
title('Permeability reduction (after phase II)', 'FontSize', fS, ...
                              'FontName', 'Arial', 'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', -aquiferLength : 25 : aquiferLength, ...
    'YTick', -aquiferWidth : 10 : aquiferWidth, 'color', 'none', ...
    'FontName', 'Arial');
ylim([-aquiferWidth, aquiferWidth]);
rectangle('Position', [10, -aquiferWidth, 5, 2 * aquiferWidth], ...
              'LineWidth', lW, 'LineStyle', '-', 'edgecolor', '[0 0 0]');

% print -depsc2 Fig7.eps
