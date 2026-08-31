%% Three-dimensional MICP treatment of a leakage system
%
% Set up and solve the three-dimensional flow and leakage system (3Dfls).
% The example first simulates MICP treatment and then assesses CO2 leakage
% before treatment and after three treatment stages.
%
% In MATLAB, this example produces Figures 12 and 13 from:
%
% Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E., 2021. Practical
% approaches to study microbially induced calcite precipitation at the
% field scale. International Journal of Greenhouse Gas Control 106, 103256.
%
% https://doi.org/10.1016/j.ijggc.2021.103256
%
% In GNU Octave, the example writes the results to the
% `vtk_micp_3Dfls` directory for visualization in ParaView.
%
% The example assumes that MRST, DistMesh, and the `ad-micp` module are
% available on the MATLAB or GNU Octave path.

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

L = 500;  % Reservoir half-length and half-width, m
H = 160;  % Reservoir height, m

if exist('AD_MICP_TEST', 'var')
    % Use fewer vertical layers in tests to reduce memory consumption.
    nz = 0.1;
else
    nz = 0.5;
end

% Fix the random seed to improve mesh reproducibility.
rand('seed', 0);
randn('seed', 0);

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    if exist('AD_MICP_TEST', 'var')
        % Use a coarser grid to reduce memory consumption.
        minimumGridSize = 0.6;
        mediumGridSize = 30;
        maximumGridSize = 500;
    else
        minimumGridSize = 0.3;
        mediumGridSize = 15;
        maximumGridSize = 500;
    end

    distanceFunction = @(points) ...
        drectangle(points, -L, L, -L, L);

    sizeFunction = @(points) min( ...
        min(1 + 0.58 .* abs(dcircle(points, -100, 0, 0)), ...
        mediumGridSize) .* ...
        (abs(dcircle(points, -100, 0, 0)) < 50) + ...
        min(mediumGridSize + ...
        0.6 .* abs(dcircle(points, -100, 0, 50)), ...
        maximumGridSize) .* ...
        (abs(dcircle(points, -100, 0, 0)) >= 50), ...
        min(minimumGridSize + ...
        0.6 .* abs(dcircle(points, 0, 0, 0)), ...
        mediumGridSize) .* ...
        (abs(dcircle(points, 0, 0, 0)) < 50) + ...
        min(mediumGridSize + ...
        0.6 .* abs(dcircle(points, 0, 0, 50)), ...
        maximumGridSize) .* ...
        (abs(dcircle(points, 0, 0, 0)) >= 50));

    boundingBox = [-L, -L; L, L];

    fixedPoints = [ ...
        -L, -L; ...
         L, -L; ...
        -L,  L; ...
         L,  L; ...
         0,  0];
    
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

    G = makeLayeredGrid( ...
        pebi(triangleGrid(gridPoints, gridTriangles)), nz .* H);
else
    maximumGridSize = 30;

    distanceFunction = @(points) ...
        drectangle(points, -L, L, -L, L);

    boundingBox = [-L, -L; L, L];

    fixedPoints = [ ...
        -L, -L; ...
         L, -L; ...
        -L,  L; ...
         L,  L];

    [gridPoints, ~] = distmesh2d( ...
        distanceFunction, @huniform, maximumGridSize, ...
        boundingBox, fixedPoints);

    close

    wellPoints = [];

    for radius = 50 .* exp(-3 : 0.125 : 0)
        [xCoordinates, yCoordinates, ~] = cylinder(radius, 28);
        wellPoints = [wellPoints, ...
                      [xCoordinates(1, :); yCoordinates(1, :)]]; %#ok<AGROW>
    end

    wellPoints = [wellPoints, [0; 0]];
    injectionWellPoints = bsxfun(@plus, wellPoints, [-100; 0]);

    leakagePoints = [];

    for radius = 50 .* exp(-5.1 : 0.125 : 0)
        [xCoordinates, yCoordinates, ~] = cylinder(radius, 28);
        leakagePoints = [leakagePoints, ...
                         [xCoordinates(1, :); yCoordinates(1, :)]]; %#ok<AGROW>
    end

    leakagePoints = [leakagePoints, [0; 0]];

    gridPoints = unique([ ...
        injectionWellPoints'; ...
        leakagePoints'; ...
        gridPoints(:, 1), gridPoints(:, 2); ...
        0, 0], ...
        'rows');

    G = makeLayeredGrid( ...
        pebi(triangleGrid(gridPoints)), nz .* H);
end

%% Adjust vertical grid coordinates

lowerExponentStart = -4;
lowerExponentStep = 4 ./ (nz .* 30);
lowerExponents = lowerExponentStart : lowerExponentStep : 0;

upperExponentStart = -4;
upperExponentStep = 4 ./ (nz .* 30);
upperExponents = upperExponentStart : upperExponentStep : 0;

numberOfNodesPerLayer = ...
    G.nodes.num ./ (nz .* H + 1);

lowerCoordinates = 30 .* exp(upperExponents);

for layerIndex = 0 : nz .* 30
    firstNode = 1 + numberOfNodesPerLayer .* layerIndex;
    lastNode = numberOfNodesPerLayer .* (layerIndex + 1);

    G.nodes.coords(firstNode : lastNode, 3) = ...
        ones(numberOfNodesPerLayer, 1) .* ...
        lowerCoordinates(layerIndex + 1);
end

for layerIndex = nz .* 30 + 1 : nz .* 130 - 1
    firstNode = 1 + numberOfNodesPerLayer .* layerIndex;
    lastNode = numberOfNodesPerLayer .* (layerIndex + 1);

    G.nodes.coords(firstNode : lastNode, 3) = ...
        G.nodes.coords(firstNode : lastNode, 3) ./ nz;
end

for layerIndex = nz .* 130 : nz .* 160
    firstNode = 1 + numberOfNodesPerLayer .* layerIndex;
    lastNode = numberOfNodesPerLayer .* (layerIndex + 1);

    G.nodes.coords(firstNode : lastNode, 3) = ...
        (G.nodes.coords(firstNode : lastNode, 3) - nz .* 130) .* ...
        exp(lowerExponents(1 + layerIndex - nz .* 130)) ./ nz + 130;
end

G = computeGeometry(G);
cellCentroids = G.cells.centroids;

%% Remove caprock cells outside the leakage path

leakageHalfWidth = 50 .* exp(-5.1) ./ 2;

cellsToRemove = ...
    (cellCentroids(:, 1) < -leakageHalfWidth | ...
    cellCentroids(:, 1) > leakageHalfWidth) & ...
    (cellCentroids(:, 3) < 130 & ...
    cellCentroids(:, 3) > 30);

G = removeCells(G, cellsToRemove);

cellCentroids = G.cells.centroids;

cellsToRemove = ...
    (cellCentroids(:, 3) < 130 & ...
    cellCentroids(:, 3) > 30) & ...
    (cellCentroids(:, 2) < -0.1 | ...
    cellCentroids(:, 2) > 0.1);

G = removeCells(G, cellsToRemove);
G = computeGeometry(G);

cellCentroids = G.cells.centroids;
cellTemplate = ones(G.cells.num, 1);

% Preserve the original variables used by the plotting sections.
c = cellCentroids;
C = cellTemplate;

%% Rock properties

K0 = 2e-14 .* cellTemplate; % Aquifer permeability, m^2

activeCellMap = G.cells.indexMap;

leakageCells = activeCellMap( ...
    cellCentroids(:, 1) > -leakageHalfWidth & ...
    cellCentroids(:, 1) < leakageHalfWidth & ...
    cellCentroids(:, 2) < leakageHalfWidth & ...
    cellCentroids(:, 2) > -leakageHalfWidth);

isLeakageCell = ismember(activeCellMap, leakageCells);

K0(isLeakageCell) = 1e-12; % Leakage-path permeability, m^2

porosity = 0.15;           % Aquifer and leakage porosity, [-]

rock = makeRock(G, K0, porosity);

%% Fluid and MICP model properties

fluid.muw = 2.535e-4;       % Water viscosity, Pa s
fluid.muO = 3.95e-5;        % CO2 viscosity, Pa s
fluid.bW = @(pressure) ...
    0 .* pressure + 1;      % Water formation volume factor, [-]
fluid.bO = @(pressure) ...
    0 .* pressure + 1;      % CO2 formation volume factor, [-]
fluid.rhoWS = 1045;         % Water density, kg/m^3
fluid.rhoOS = 479;          % CO2 density, kg/m^3

fluid.rho_b = 35;           % Biofilm density, kg/m^3
fluid.rho_c = 2710;         % Calcite density, kg/m^3
fluid.k_str = 2.6e-10;      % Detachment coefficient, m/(Pa s)
fluid.diffm = 0;            % Microorganism diffusion, m^2/s
fluid.diffo = 0;            % Oxygen diffusion, m^2/s
fluid.diffu = 0;            % Urea diffusion, m^2/s
fluid.alphaL = 0;           % Longitudinal dispersivity, m
fluid.alphaT = 0;           % Transverse dispersivity, m
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
    (porosity - fluid.crit), 0);

fluid.K = @(currentPorosity) max( ...
    (K0 .* normalizedPorosity(currentPorosity) .^ fluid.eta + ...
    fluid.kmin) .* K0 ./ (K0 + fluid.kmin), fluid.kmin);

%% MICP injection strategy
%
% Each row of `M` contains:
%
%   1. Phase duration, s
%   2. Simulation timestep, s
%   3. Injection rate, m^3/s
%   4. Microorganism concentration, kg/m^3
%   5. Oxygen concentration, kg/m^3
%   6. Urea concentration, kg/m^3

dt = hour;

M = [ ...
     15 * hour, dt, 3, 0.01, 0,    0; ...
     11 * hour, dt, 3, 0,    0,    0; ...
     74 * hour, dt, 0, 0,    0,    0; ...
     30 * hour, dt, 3, 0,    0.04, 0; ...
      5 * hour, dt, 3, 0,    0,    0; ...
     25 * hour, dt, 0, 0,    0,    0; ...
     40 * hour, dt, 3, 0,    0,  300; ...
     10 * hour, dt, 3, 0,    0,    0; ...
    390 * hour, dt, 0, 0,    0,    0; ...
     30 * hour, dt, 3, 0,    0.04, 0; ...
     20 * hour, dt, 3, 0,    0,    0; ...
     20 * hour, dt, 0, 0,    0,    0; ...
     20 * hour, dt, 3, 0,    0,  300; ...
     20 * hour, dt, 3, 0,    0,    0; ...
     90 * hour, dt, 0, 0,    0,    0; ...
     20 * hour, dt, 3, 0,    0,  300; ...
     20 * hour, dt, 3, 0,    0,    0; ...
    110 * hour, dt, 0, 0,    0,    0];

if exist('AD_MICP_TEST', 'var')
    dt = 1;

    M(M(:, 3) > 0, 3) = 1e-1;
    M(:, 1) = dt;
    M(:, 2) = dt;

    M(10, 5) = 0.001;
    M(7, 6) = 1;
    M(13, 6) = 1;
    M(16, 6) = 1;
end

N = size(M, 1);

%% MICP injection well

r = 0.15;
Whu = 1 / 10;
Whb = 1 - Whu;

[~, iw] = min( ...
    (cellCentroids(:, 1) + 100) .^ 2 + ...
    cellCentroids(:, 2) .^ 2);

allCellIndices = (1 : G.cells.num)';

cellsWu = allCellIndices( ...
    abs(cellCentroids(:, 1) - cellCentroids(iw, 1)) < 0.01 & ...
    abs(cellCentroids(:, 2) - cellCentroids(iw, 2)) < 0.01 & ...
    cellCentroids(:, 3) > 130 & ...
    cellCentroids(:, 3) < 133);

W = addWell([], G, rock, cellsWu, 'Type', 'rate', ...
    'Comp_i', [1, 0], 'Val', Whu .* M(1, 3), 'Radius', r);

cellsWb = allCellIndices( ...
    abs(cellCentroids(:, 1) - cellCentroids(iw, 1)) < 0.01 & ...
    abs(cellCentroids(:, 2) - cellCentroids(iw, 2)) < 0.01 & ...
    cellCentroids(:, 3) > 133);

W = addWell(W, G, rock, cellsWb, 'Type', 'rate', ...
    'Comp_i', [1, 0], 'Val', Whb .* M(1, 3), 'Radius', r);

for wellIndex = 1 : numel(W)
    W(wellIndex).m = 0;
    W(wellIndex).o = 0;
    W(wellIndex).u = 0;
end

W(1).m = M(1, 4);
W(1).o = M(1, 5);
W(1).u = M(1, 6);

G.injectionwellonboundary = 0;

%% Gravity

gravity reset
gravity on

gravityMagnitude = norm(gravity);

%% Boundary conditions

boundaryFaceIndices = boundaryFaces(G);

openBoundaryFaces = boundaryFaceIndices( ...
    abs(G.faces.normals(boundaryFaceIndices, 1)) > eps & ...
    (G.faces.centroids(boundaryFaceIndices, 1) < -L + 2 | ...
    G.faces.centroids(boundaryFaceIndices, 1) > L - 2));

boundaryPressure = ...
    G.faces.centroids(openBoundaryFaces, 3) .* ...
    fluid.rhoWS .* gravityMagnitude;

bc = addBC([], openBoundaryFaces, 'pressure', boundaryPressure, ...
                                           'sat', [0 0]);

numberOfBoundaryFaces = size(bc.sat, 1);

bc.m = zeros(numberOfBoundaryFaces, 1);
bc.o = zeros(numberOfBoundaryFaces, 1);
bc.u = zeros(numberOfBoundaryFaces, 1);
bc.b = zeros(numberOfBoundaryFaces, 1);
bc.c = zeros(numberOfBoundaryFaces, 1);

%% Construct MICP schedule

stepsPerPhase = round(M(:, 1) ./ M(:, 2));
nt = sum(stepsPerPhase);

timesteps = zeros(nt, 1);
controlIndices = zeros(nt, 1);

firstStep = 1;

for phaseIndex = 1 : N
    lastStep = firstStep + stepsPerPhase(phaseIndex) - 1;

    timesteps(firstStep : lastStep) = M(phaseIndex, 2);
    controlIndices(firstStep : lastStep) = phaseIndex;

    firstStep = lastStep + 1;
end

schedule = simpleSchedule(timesteps, 'W', W, 'bc', bc);

baseControl = schedule.control(1);
schedule.control = repmat(baseControl, N, 1);

for phaseIndex = 1 : N
    schedule.control(phaseIndex).W(1).val = Whu .* M(phaseIndex, 3);
    schedule.control(phaseIndex).W(2).val = Whb .* M(phaseIndex, 3);
    schedule.control(phaseIndex).W(1).m = M(phaseIndex, 4);
    schedule.control(phaseIndex).W(1).o = M(phaseIndex, 5);
    schedule.control(phaseIndex).W(1).u = M(phaseIndex, 6);
end

schedule.step.control = controlIndices;
schedule.step.val = timesteps;

phaseEndSteps = cumsum(stepsPerPhase);

% Maximum injected oxygen and urea concentrations
fluid.Comax = max(M(:, 5));
fluid.Cumax = max(M(:, 6));

%% MICP model and initial state

model = MICPModel(G, rock, fluid);

initialPressure = ...
    cellCentroids(:, 3) .* fluid.rhoWS .* gravityMagnitude;

state0 = initState(G, W, initialPressure, [1, 0]);

state0.m = zeros(G.cells.num, 1);
state0.o = zeros(G.cells.num, 1);
state0.u = zeros(G.cells.num, 1);
state0.b = zeros(G.cells.num, 1);
state0.c = zeros(G.cells.num, 1);

%% Simulate MICP treatment

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    [~, states] = simulateScheduleAD(state0, model, schedule);
else
    micpPlotFunction = getPlotAfterStepMICP(state0, model, 340, 20);

    [~, states] = simulateScheduleAD(state0, model, schedule, ...
                                      'afterStepFn', micpPlotFunction);
end

%% Select treatment states for CO2 assessment

statesa = state0;

phaseIIndex = min(9, N);
phaseIIIndex = min(15, N);

statesb = states{phaseEndSteps(phaseIIndex)};
statesc = states{phaseEndSteps(phaseIIIndex)};
statesd = states{end};

%% CO2 assessment schedule

co2TimeStep = day;
ntco2 = 100;

if exist('AD_MICP_TEST', 'var')
    co2TimeStep = minute;
    ntco2 = 2;
end

co2Timesteps = repmat(co2TimeStep, ntco2, 1);

%% CO2 injection well

QCO2 = 1600 ./ day;

co2WellCells = allCellIndices( ...
    abs(cellCentroids(:, 1) - cellCentroids(iw, 1)) < 0.01 & ...
    abs(cellCentroids(:, 2) - cellCentroids(iw, 2)) < 0.01 & ...
    cellCentroids(:, 3) > 130);

co2Well = addWell([], G, rock, co2WellCells, 'Type', 'rate', ...
    'Comp_i', [0, 1], 'Val', QCO2, 'Radius', r);

co2Schedule = simpleSchedule( ...
    co2Timesteps, 'W', co2Well, 'bc', bc);

co2InitialState = initState( ...
    G, co2Well, initialPressure, [1, 0]);

%% Build rock states before and after MICP treatment

porositya = max(porosity - statesa.c - statesa.b, ...
                model.minimumPorosity);
porosityf = max(porosity - statesb.c - statesb.b, ...
                model.minimumPorosity);
porosityg = max(porosity - statesc.c - statesc.b, ...
                model.minimumPorosity);
porosityh = max(porosity - statesd.c - statesd.b, ...
                model.minimumPorosity);

rocka = makeRock(G, max(fluid.K(porositya), fluid.kmin), porositya);
rockb = makeRock(G, max(fluid.K(porosityf), fluid.kmin), porosityf);
rockc = makeRock(G, max(fluid.K(porosityg), fluid.kmin), porosityg);
rockd = makeRock(G, max(fluid.K(porosityh), fluid.kmin), porosityh);

modela = CO2Model(G, rocka, fluid);
modelb = CO2Model(G, rockb, fluid);
modelc = CO2Model(G, rockc, fluid);
modeld = CO2Model(G, rockd, fluid);

%% Identify leakage face

allFaceIndices = (1 : G.faces.num)';

leakageFaces = allFaceIndices( ...
    G.faces.centroids(:, 3) < 80 + 1000 .* eps & ...
    G.faces.centroids(:, 3) > 80 - 1000 .* eps);

%% Simulate CO2 migration

lr0 = zeros(ntco2, 1);
lr1 = zeros(ntco2, 1);
lr2 = zeros(ntco2, 1);
lr3 = zeros(ntco2, 1);

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    [~, statesco2] = simulateScheduleAD( ...
        co2InitialState, modela, co2Schedule);

    for stepIndex = 1 : ntco2
        lr0(stepIndex) = abs( ...
            statesco2{stepIndex}.flux(leakageFaces(1), 2));
    end

    statese = statesco2{end};
    clear statesco2

    [~, statesco2] = simulateScheduleAD( ...
        co2InitialState, modelb, co2Schedule);

    for stepIndex = 1 : ntco2
        lr1(stepIndex) = abs( ...
            statesco2{stepIndex}.flux(leakageFaces(1), 2));
    end

    statesf = statesco2{end};
    clear statesco2

    [~, statesco2] = simulateScheduleAD( ...
        co2InitialState, modelc, co2Schedule);

    for stepIndex = 1 : ntco2
        lr2(stepIndex) = abs( ...
            statesco2{stepIndex}.flux(leakageFaces(1), 2));
    end

    statesg = statesco2{end};
    clear statesco2

    [~, statesco2] = simulateScheduleAD( ...
        co2InitialState, modeld, co2Schedule);

    for stepIndex = 1 : ntco2
        lr3(stepIndex) = abs( ...
            statesco2{stepIndex}.flux(leakageFaces(1), 2));
    end

    statesh = statesco2{end};
    clear statesco2
else
    co2PlotFunctionA = getPlotAfterStepCO2( ...
        co2InitialState, modela, 340, 20);
    co2PlotFunctionB = getPlotAfterStepCO2( ...
        co2InitialState, modelb, 340, 20);
    co2PlotFunctionC = getPlotAfterStepCO2( ...
        co2InitialState, modelc, 340, 20);
    co2PlotFunctionD = getPlotAfterStepCO2( ...
        co2InitialState, modeld, 340, 20);

    [~, statesco2] = simulateScheduleAD( ...
        co2InitialState, modela, co2Schedule, ...
        'afterStepFn', co2PlotFunctionA);

    for stepIndex = 1 : ntco2
        lr0(stepIndex) = abs( ...
            statesco2{stepIndex}.flux(leakageFaces(1), 2));
    end

    statese = statesco2{end};
    clear statesco2

    [~, statesco2] = simulateScheduleAD( ...
        co2InitialState, modelb, co2Schedule, ...
        'afterStepFn', co2PlotFunctionB);

    for stepIndex = 1 : ntco2
        lr1(stepIndex) = abs( ...
            statesco2{stepIndex}.flux(leakageFaces(1), 2));
    end

    statesf = statesco2{end};
    clear statesco2

    [~, statesco2] = simulateScheduleAD( ...
        co2InitialState, modelc, co2Schedule, ...
        'afterStepFn', co2PlotFunctionC);

    for stepIndex = 1 : ntco2
        lr2(stepIndex) = abs( ...
            statesco2{stepIndex}.flux(leakageFaces(1), 2));
    end

    statesg = statesco2{end};
    clear statesco2

    [~, statesco2] = simulateScheduleAD( ...
        co2InitialState, modeld, co2Schedule, ...
        'afterStepFn', co2PlotFunctionD);

    for stepIndex = 1 : ntco2
        lr3(stepIndex) = abs( ...
            statesco2{stepIndex}.flux(leakageFaces(1), 2));
    end

    statesh = statesco2{end};
    clear statesco2
end

% Preserve the original timestep variable used by the Figure 13b plotting
% section.
dt = co2TimeStep;

%% Export results in GNU Octave

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    outputDirectory = fullfile(pwd, 'vtk_micp_3Dfls');

    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    micpOutput = fullfile(outputDirectory, 'states_3Dfls');
    co2Output = fullfile(outputDirectory, 'states_3Dfls_CO2');

    mrsttovtk(G, states, micpOutput, '%f');
    mrsttovtk(G, statesh, co2Output, '%f');

    fprintf(['VTK results written to:\n' ...
             '  %s.pvd\n' ...
             '  %s.pvd\n'], ...
             micpOutput, co2Output);

    return
end
 
% Figure 12 paper (MATLAB)
porosityf = porosity - statesb.c - statesb.b;
porosityg = porosity - statesc.c - statesc.b;
porosityh = porosity - statesd.c - statesd.b;

cellsF =  1 : G.cells.num;
cellsf =  1 : G.cells.num;
cellsf = cellsf(c(:, 1) < 0 & c(:, 2) < 0);
idx = ismember(cellsF, cellsf);
cellsFa =  1 : G.cells.num;
cellsfa =  1 : G.cells.num;
cellsfa = cellsfa((c(:, 1) < 0 | c(:, 2) < 0) & c(:, 3) < 30);
idxa = ismember(cellsFa, cellsfa);

c = flipud(jet);
sz = size(c, 1);
ccc = c((round(70 * sz / 256)) : end, :);
c = c((round(70 * sz / 256)) : (round(100 * sz / 256)), :);
cc(:, 1) = (0.75 : 0.01 : 1)';
cc(:, 2) = (0.75 : 0.01 : 1)';
cc(:, 3) = (0.75 : -0.03 : 0)';
fS = 8;
lW = 1;
figure;
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 9.11 4.83]);
hold on
n1 = subplot(2, 4, 1);
view(340, 45);
xlim([-L L])
ylim([-L L])
zlim([0 H])
xlabel({'x [m]' ; '(a)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('y [m]', 'FontSize', fS, 'FontName', 'Arial');
zlabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
plotGrid(G, ~idxa' .* K0 > 0, 'FaceColor', '[0.75 0.75 0.75]');
s = plotCellData(G, K0, ~idxa' .* K0 > 0);
title('Initial permeability', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', -L : 250 : L, 'YTick', -L : 250 : L, ...
             'ZTick', -L : 40 : 160, 'color', 'none', 'FontName', 'Arial');
colormap (n1, cc);
caxis([2e-14 1e-12]);
cb = colorbar; 
title(cb, 'm$^2$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.26 0.71 0.005 0.08], 'YTick', [2e-14 1e-12]);
line([-175 0], [0 0], [-10 30], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([-175 0], [0 0], [30 130], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
ax1=axes('position', [0.15 0.81 0.03 0.03], 'YAxisLocation', 'right');
box on
axis([-0.25 0.25 0 160]);
xlim([-0.2 0.2])
zlim([30 130])
s = plotCellData(G, K0);
s.EdgeColor = 'none';
colormap (ax1, cc);
caxis([2e-14 1e-12]);
set(gca, 'FontSize', 6, 'XTick', [-0.15 0.15], 'ZTick', [30 130], ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 0)
n2 = subplot(2, 4, 2);
view(340, 20);
colormap (n2, ccc);
caxis([0 100]);
cb = colorbar; 
title(cb, '$\%$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.47 0.675 0.005 0.08], 'YTick', [0 50 100]);
xlim([-L L])
ylim([-L L])
zlim([0 H])
xlabel({'x [m]' ; '(b)'}, 'FontSize', fS, 'FontName', 'Arial');
zlabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
ylabel('y [m]', 'FontSize', fS, 'FontName', 'Arial');
plotGrid(G, idx, 'FaceColor', 'none', 'EdgeAlpha', 0.25);
s = plotCellData(G, ~idx' .* 100 .* (K0 - fluid.K(porosity - statesb.c ...
                                     - statesb.b)) ./ K0, ~idx' .* K0 > 0);
s.EdgeColor = 'none';
title('Permeability (phase I MICP)', 'FontSize', fS, ...
                              'FontName', 'Arial', 'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', -L : 250 : L, 'YTick', -L : 250 : L, ...
              'ZTick', 0 : 40 : 160, 'color', 'none', 'FontName', 'Arial');
view(340, 20);
set(gca, 'FontName', 'Arial');
line([-175 0], [0 0], [65 30], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([-175 0], [0 0], [90 130], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
ax2 = axes('position', [0.355 0.74 0.03 0.03], 'YAxisLocation', 'left');
box on
axis([-0.25 0.25 0 160]);
xlim([-0.2 0.2])
zlim([30 130])
s = plotCellData(G, ~idx' .* 100 .* (K0 - fluid.K(porosity - statesb.c ...
                                     - statesb.b)) ./ K0, ~idx' .* K0 > 0);
s.EdgeColor = 'none';
set(gca, 'FontSize', 6, 'XTick', ([-0.15 0.15]), 'ZTick', [30 130], ...
                                     'color', 'none', 'FontName', 'Arial');
view(0,0)
colormap (ax2, ccc);
caxis([0 100]);

n3 = subplot(2, 4, 3);
view(340, 20);
colormap (n3, ccc);
caxis([0 100]);
cb = colorbar; 
title(cb, '$\%$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.677 0.675 0.005 0.08], 'YTick', [0 50 100]);
xlabel({'x [m]' ; '(c)'}, 'FontSize', fS, 'FontName', 'Arial');
zlabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
ylabel('y [m]', 'FontSize', fS, 'FontName', 'Arial');
plotGrid(G, idx, 'FaceColor', 'none', 'EdgeAlpha', 0.25);
s = plotCellData(G, ~idx' .* 100 .* (K0 - fluid.K(porosity - statesc.c ...
                                     - statesc.b)) ./ K0, ~idx' .* K0 > 0);
s.EdgeColor = 'none';
title('Permeability (phase II MICP)', 'FontSize', fS, ...
                              'FontName', 'Arial', 'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', -L : 250 : L, 'YTick', -L : 250 : L, ...
              'ZTick', 0 : 40 : 160, 'color', 'none', 'FontName', 'Arial');
view(340, 20);
set(gca, 'FontName', 'Arial');
line([-175 0], [0 0], [65 30], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([-175 0], [0 0], [90 130], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
ax3 = axes('position', [0.561 0.74 0.03 0.03], 'YAxisLocation', 'left');
box on
axis([-0.25 0.25 0 160]);
xlim([-0.2 0.2])
zlim([30 130])
s = plotCellData(G, ~idx' .* 100 .* (K0 - fluid.K(porosity - statesc.c ...
                                     - statesc.b)) ./ K0, ~idx' .* K0 > 0);
s.EdgeColor = 'none';
set(gca, 'FontSize', 6, 'XTick', ([-0.15 0.15]), 'ZTick', [30 130], ...
                                     'color', 'none', 'FontName', 'Arial');
view(0,0)
colormap (ax3, ccc);
caxis([0 100]);

n4 = subplot(2, 4, 4);
view(340, 20);
colormap (n4, ccc);
caxis([0 100]);
cb = colorbar; 
title(cb, '$\%$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.88 0.675 0.005 0.08], 'YTick', [0 50 100]);
xlabel({'x [m]' ; '(d)'}, 'FontSize', fS, 'FontName', 'Arial');
zlabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
ylabel('y [m]', 'FontSize', fS, 'FontName', 'Arial');
plotGrid(G, idx, 'FaceColor', 'none', 'EdgeAlpha', 0.25);
s = plotCellData(G, ~idx' .* 100 .* (K0 - fluid.K(porosity - statesd.c ...
                                   - statesd.b)) ./ K0, ~idx' .* K0 > 0);
s.EdgeColor = 'none';
title('Permeability (phase II MICP)', 'FontSize', fS, ...
                              'FontName', 'Arial', 'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', -L : 250 : L, 'YTick', -L : 250 : L, ...
              'ZTick', 0 : 40 : 160, 'color', 'none', 'FontName', 'Arial');
view(340, 20);
set(gca, 'FontName', 'Arial');
line([-175 0], [0 0], [65 30], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([-175 0], [0 0], [90 130],'Color','black','LineStyle','--', ...
                                                          'LineWidth', lW);
ax3 = axes('position', [0.766 0.74 0.03 0.03], 'YAxisLocation', 'left');
box on
axis([-0.25 0.25 0 160]);
xlim([-0.2 0.2])
zlim([30 130])
s = plotCellData(G, ~idx' .* 100 .* (K0 - fluid.K(porosity - statesd.c ...
                                     - statesd.b)) ./ K0, ~idx' .* K0 > 0);
s.EdgeColor = 'none';
set(gca, 'FontSize', 6, 'XTick', ([-0.15 0.15]), 'ZTick', [30 130], ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 0)
colormap (ax3, ccc);
caxis([0 100]);

n5 = subplot(2, 4, 5);
view(340, 20);
colormap (n5, c);
caxis([0 75]);
cb = colorbar; 
title(cb, 'kg/m$^3$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.26 0.2 0.005 0.08], 'YTick', [0 25 50 75]);
xlabel({'x [m]' ; '(e)'}, 'FontSize', fS, 'FontName', 'Arial');
zlabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
ylabel('y [m]', 'FontSize', fS, 'FontName', 'Arial');
plotGrid(G, idx, 'FaceColor', 'none', 'EdgeAlpha', 0.25);
s = plotCellData(G, ~idx' .* porosity .* fluid.rhoOS .* ...
                                         statese.s(:, 2), ~idx' .* K0 > 0);
s.EdgeColor = 'none';
title('CO$_2$ (100 days)', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', -L : 250 : L, 'YTick', -L : 250 : L, ...
              'ZTick', 0 : 40 : 160, 'color', 'none', 'FontName', 'Arial');
set(gca, 'FontName', 'Arial');
line([-175 0], [0 0], [65 30], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([-175 0], [0 0], [90 130], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
ax4 = axes('position', [0.15 0.27 0.03 0.03], 'YAxisLocation', 'right');
box on
axis([-0.25 0.25 0 160]);
xlim([-0.2 0.2])
zlim([30 130])
s = plotCellData(G, ~idx' .* porosity .* fluid.rhoOS .* ...
                                          statese.s(:, 2),~idx' .* K0 > 0);
s.EdgeColor = 'none';
colormap (ax4, c);
caxis([0 75]);
set(gca, 'FontSize', 6, 'XTick', [-0.15 0.15], 'ZTick', [30 130], ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 0)
n6 = subplot(2, 4, 6);
view(340, 20);
colormap (n6, c);
caxis([0 75]);
cb = colorbar; 
title(cb, 'kg/m$^3$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.47 0.2 0.005 0.08], 'YTick', [0 25 50 75]);
xlabel({'x [m]' ; '(f)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('y [m]', 'FontSize', fS, 'FontName', 'Arial');
zlabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
plotGrid(G, idx, 'FaceColor', 'none', 'EdgeAlpha', 0.25);
s = plotCellData(G, ~idx' .* fluid.rhoOS .* porosityf.* ...
                                          statesf.s(:, 2),~idx' .* K0 > 0);
s.EdgeColor = 'none';
title('CO$_2$ (phase I MICP)', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', 8, 'XTick', -L : 250 : L, 'YTick', -L : 250 : L, ...
              'ZTick', 0 : 40 : 160, 'color', 'none', 'FontName', 'Arial');
set(gca, 'FontName', 'Arial');
line([-175 0], [0 0], [65 30], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([-175 0], [0 0], [90 130], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
ax5 = axes('position', [0.355 0.27 0.03 0.03], 'YAxisLocation', 'right');
box on
axis([-0.25 0.25 0 160]);
xlim([-0.2 0.2])
zlim([30 130])
s = plotCellData(G, fluid.rhoOS .* porosityf .* statesf.s(:, 2));
s.EdgeColor = 'none';
colormap (ax5, c);
caxis([0 75]);
set(gca, 'FontSize', 6, 'XTick', [-0.15 0.15], 'ZTick', [30 130], ...
                                     'color', 'none', 'FontName', 'Arial');
view(0,0)

n7 = subplot(2, 4, 7);
view(340, 20);
colormap (n7, c);
caxis([0 75]);
cb = colorbar; 
title(cb, 'kg/m$^3$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.677 0.2 0.005 0.08], 'YTick', [0 25 50 75]);
xlabel({'x [m]' ; '(g)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('y [m]', 'FontSize', fS, 'FontName', 'Arial');
zlabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
plotGrid(G, idx, 'FaceColor', 'none', 'EdgeAlpha', 0.25);
s = plotCellData(G, ~idx' .* fluid.rhoOS .* porosityg .* ...
                                          statesg.s(:, 2),~idx' .* K0 > 0);
s.EdgeColor = 'none';
title('CO$_2$ (phase II MICP)', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', -L : 250 : L, 'YTick', -L : 250 : L, ...
              'ZTick', 0 : 40 : 160, 'color', 'none', 'FontName', 'Arial');
set(gca, 'FontName', 'Arial');
line([-175 0], [0 0], [65 30], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([-175 0], [0 0], [90 130], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
ax6 = axes('position', [0.561 0.27 0.03 0.03], 'YAxisLocation', 'right');
box on
axis([-0.25 0.25 0 160]);
xlim([-0.2 0.2])
zlim([30 130])
s = plotCellData(G, fluid.rhoOS .* porosityg .* statesg.s(:, 2));
s.EdgeColor = 'none';
colormap (ax6, c);
caxis([0 75]);
set(gca, 'FontSize', 6,'XTick', [-0.15 0.15], 'ZTick', [30 130], ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 0)
n8 = subplot(2, 4, 8);
view(340, 20);
colormap (n8, c);
caxis([0 75]);
cb = colorbar; 
title(cb, 'kg/m$^3$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.88 0.2 0.005 0.08], 'YTick', [0 25 50 75]);
xlabel({'x [m]' ; '(h)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('y [m]', 'FontSize', fS, 'FontName', 'Arial');
zlabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
plotGrid(G, idx, 'FaceColor', 'none', 'EdgeAlpha', 0.25);
s = plotCellData(G, ~idx' .* fluid.rhoOS .* porosityh .* ...
                                         statesh.s(:, 2), ~idx' .* K0 > 0);
s.EdgeColor = 'none';
title('CO$_2$ (phase III MICP)', 'FontSize', fS, 'FontName', 'Arial', ...
                                                    'Interpreter','latex');
set(gca, 'FontSize', fS, 'XTick', -L : 250 : L, 'YTick', -L : 250 : L, ...
              'ZTick', 0 : 40 : 160, 'color', 'none', 'FontName', 'Arial');
set(gca, 'FontName', 'Arial');
line([-175 0], [0 0], [65 30], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([-175 0], [0 0], [90 130], 'Color', 'black', 'LineStyle','--', ...
                                                          'LineWidth', lW);
ax6 = axes('position', [0.766 0.27 0.03 0.03], 'YAxisLocation', 'right');
box on
axis([-0.25 0.25 0 160]);
xlim([-0.2 0.2])
zlim([30 130])
s = plotCellData(G, fluid.rhoOS .* porosityh .* statesh.s(:, 2));
s.EdgeColor = 'none';
colormap (ax6, c);
caxis([0 75]);
set(gca, 'FontSize', 6, 'XTick', [-0.15 0.15], 'ZTick', [30 130], ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 0)
%print -depsc2 Fig12.eps

% Figures 13a and 13b paper
clear m o u b c K vc v cell_leak
cells = 1 : G.cells.num;
nt_micp = cumsum(schedule.step.val) / hour;
cell_leak = cells(G.cells.centroids(:, 3) < 130 & ...
                                             G.cells.centroids(:, 3) > 30);

for i = 1 : nt  
    m(i) = mean(states{i}.m(cell_leak));
    o(i) = mean(states{i}.o(cell_leak));
    u(i) = mean(states{i}.u(cell_leak));
    b(i) = mean(states{i}.b(cell_leak));
    c(i) = mean(states{i}.c(cell_leak));
    currentPorosity = max(porosity - states{i}.c - states{i}.b, ...
                      model.minimumPorosity);
    Ki = fluid.K(currentPorosity);
    K(i) = mean(Ki(cell_leak) ./ K0(cell_leak));
    vc = faceFlux2cellVelocity(G, states{i}.flux(:));
    v(i) = mean(sqrt(sum(vc(cell_leak, :) .^ 2, 2)));
end
lW = 3; 
fS = 11;
figure('Units', 'inches', 'PaperPosition', [0 0 6.83 6.83]);
hold on
plot(nt_micp, v / max(v), 'color', [0 0.74 1], 'LineWidth', lW, ...
                                                         'LineStyle', '-');
plot(nt_micp, m / max(m), 'color', [0 0.8 0], 'LineWidth', lW, ...
                                                         'LineStyle', '-');
plot(nt_micp, o / max(o), 'color', [1 0.5 0.9], 'LineWidth', lW, ...
                                                         'LineStyle', '-');
plot(nt_micp, u / max(u), 'color', [1 0.9 0], 'LineWidth', lW, ...
                                                         'LineStyle', '-');
plot(nt_micp, b / max(b), 'color', [0 0.4 0], 'LineWidth', lW, ...
                                                         'LineStyle', ':');
plot(nt_micp, c / max(c), 'color', [1 0.2 0.2], 'LineWidth', lW, ...
                                                         'LineStyle', ':');
plot(nt_micp, K, 'color', [0 0 0], 'LineWidth', lW, 'LineStyle', ':');
line([600 600], [0 1], [0 0], 'Color', [0 0 0], 'LineStyle', '--', ...
                                                           'LineWidth', 1);
line([800 800], [0 1], [0 0], 'Color', [0 0 0], 'LineStyle', '--', ...
                                                           'LineWidth', 1);                                                       
hold off

text(250, 1.02, 'Phase I', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
text(650, 1.02, 'Phase II', 'FontSize', fS,'Interpreter','latex', ...
                                                      'FontName', 'Arial');
text(850, 1.02, 'Phase III', 'FontSize', fS,'Interpreter','latex', ...
                                                      'FontName', 'Arial');                                                   
xlim([0 nt_micp(end)]);
xlabel({'Time [h]' ; '(a)'}, 'FontSize', fS, 'Interpreter', 'latex');        
ylabel('[$-$]', 'FontSize', fS, 'Interpreter', 'latex');
h = legend('$v_w/0.0153\textrm{ m/s}$', '$c_m/0.0069\textrm{ kg/m}^3$', ...
          '$c_o/0.0298\textrm{ kg/m}^3$', '$c_u/229 \textrm{ kg/m}^3$', ...
                                  '$\phi_b/0.0002$', '$\phi_c/0.0333$', ...
                                           '$K/10^{-12}\textrm{ m}^2$', ...
                                                'Interpreter', 'latex', ...
                                                           'FontSize', fS);
rect = [0.36, 0.55, 0.2, 0.25];
set(h, 'Position', rect);               
set(gca, 'FontSize', fS, 'FontName', 'Arial', 'XTick', 0 : 100 : 1000, ...
                                             'YGrid', 'on', 'XGrid', 'on');
%print -depsc2 Fig13a.eps

lW = 9;
fS = 11;
figure('Units', 'inches', 'PaperPosition', [0 0 6.83 6.83]);
hold on
plot((1 : ntco2) * dt / day, 100 * lr0 / QCO2, ...
                  'color', [1 0.2 0.2], 'LineWidth', lW, 'LineStyle', '-');
plot((1 : ntco2) * dt / day, 100 * lr1 / QCO2, ...
                     'color',[1 0.5 0], 'LineWidth', lW, 'LineStyle', '-');
plot((1 : ntco2) * dt / day, 100 * lr2 / QCO2, ...
              'color',[0.61 0.61 0.61], 'LineWidth', lW, 'LineStyle', '-');
plot((1 : ntco2) * dt / day, 100 * lr3 / QCO2, ...
                       'color',[0 0 0], 'LineWidth', lW, 'LineStyle', '-');
hold off
xlim([0 100]);
ylim([0 0.30]);
xlabel({'Time [d]' ; '(b)'}, 'FontSize', fS, 'Interpreter', 'latex');        
ylabel('CO$_2$ leakage rate/injection rate [\%]', 'FontSize', fS, ...
                                                   'Interpreter', 'latex');
grid on
legend('Without MICP', 'Phase I MICP', 'Phase II MICP', ...
                                     'Phase III MICP', 'Location', 'best');
set(gca, 'FontSize', fS, 'FontName', 'Arial', 'XTick', 0 : 20 : 100, ...
                                                 'YTick', 0 : 0.06 : 0.30);
%print -depsc2 Fig13b.eps 
