%% Two-dimensional MICP treatment of a leakage system
%
% Set up and solve the two-dimensional flow and leakage system (2Dfls).
% The example first simulates MICP treatment and then assesses CO2 leakage
% before treatment and after three treatment stages.
%
% In MATLAB, this example produces Figures 10 and 11 from:
%
% Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E., 2021. Practical
% approaches to study microbially induced calcite precipitation at the
% field scale. International Journal of Greenhouse Gas Control 106, 103256.
%
% https://doi.org/10.1016/j.ijggc.2021.103256
%
% In GNU Octave, the example writes the results to the
% `vtk_micp_2Dfls` directory for visualization in ParaView.
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

L = 500; % Reservoir length, m
H = 160; % Reservoir height, m

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    verticalCoordinates = [ ...
        0 : 0.25 : 135, ...
        H .* exp(-0.16 : 0.02 : 0)];

    horizontalCoordinates = [ ...
        0, ...
        1, ...
        50 .* exp(-3.6 : 0.25 : 0), ...
        100 - 50 .* exp(0 : -0.25 : -4.75), ...
        99.7, ...
        99.85, ...
        100.15, ...
        100.3, ...
        100.5, ...
        101, ...
        102, ...
        105 : 5 : 210, ...
        L .* exp(-0.8 : 0.05 : 0)];

    horizontalCoordinates = unique(horizontalCoordinates, 'stable');
    verticalCoordinates = unique(verticalCoordinates, 'stable');

    horizontalCoordinates(end) = L;
    verticalCoordinates(end) = H;

    G = tensorGrid( ...
        horizontalCoordinates, verticalCoordinates, [0 1]);

    G = computeGeometry(G);
    cellCentroids = G.cells.centroids;

    cellsToRemove = ...
        (cellCentroids(:, 1) < 100 - 0.3 | ...
        cellCentroids(:, 1) > 100 + 0.3) & ...
        (cellCentroids(:, 2) < 130 & ...
        cellCentroids(:, 2) > 30);

    G = removeCells(G, cellsToRemove);
    G = computeGeometry(G);
    cellCentroids = G.cells.centroids;

    cellsToRemove = ...
        (cellCentroids(:, 1) < 99.9 | ...
        cellCentroids(:, 1) > 100.1) & ...
        (cellCentroids(:, 2) < 130 & ...
        cellCentroids(:, 2) > 30);

    G = removeCells(G, cellsToRemove);
    G = computeGeometry(G);
    cellCentroids = G.cells.centroids;

    cellsToRemove = ...
        cellCentroids(:, 1) < -0.5 - eps & ...
        cellCentroids(:, 2) > 130;

    G = removeCells(G, cellsToRemove);
    G = computeGeometry(G);
    cellCentroids = G.cells.centroids;

    cellsToRemove = ...
        cellCentroids(:, 1) < 99.8 & ...
        cellCentroids(:, 2) < 30;

    G = removeCells(G, cellsToRemove);
    G = computeGeometry(G);
else
    [coarseX, coarseY] = meshgrid( ...
        [-L : 10 : -50, 180 : 10 : L], ...
        [0 : 5 : 30, 31 : 129, 130 : 5 : H]);

    [bottomX, bottomY] = meshgrid( ...
        -50 : 10 : 50, 0 : 5 : 30);

    [leakageX, leakageY] = meshgrid( ...
        [100 - 0.3, 100, 100 + 0.3], 0 : 0.25 : H);

    [leftChannelX, leftChannelY] = meshgrid( ...
        [-1, 0, 1], 130 : 0.25 : H);

    [upperLeftX, upperLeftY] = meshgrid( ...
        [-50 .* exp(0 : -0.25 : -3.6), ...
          50 .* exp(-3.6 : 0.25 : 0)], ...
        130 : 1 : H);

    [lowerLeakageX, lowerLeakageY] = meshgrid( ...
        [100 - 50 .* exp(0 : -0.25 : -5), ...
         100 + 80 .* exp(-5 : 0.125 : 0)], ...
        0 : H);

    [transitionX1, transitionY1] = meshgrid( ...
        50 .* exp(-3.6 : 0.25 : 0), 130 : 0.25 : 137.5);

    [transitionX2, transitionY2] = meshgrid( ...
        50 .* exp(-3.6 : 0.25 : 0), 137.5 : 0.5 : 145);

    [transitionX3, transitionY3] = meshgrid( ...
        100 - 50 .* exp(0 : -0.25 : -5), 130 : 0.25 : 137.5);

    [transitionX4, transitionY4] = meshgrid( ...
        100 - 50 .* exp(0 : -0.25 : -5), 137.5 : 0.5 : 145);

    gridPoints = unique([ ...
        coarseX(:), coarseY(:); ...
        bottomX(:), bottomY(:); ...
        leakageX(:), leakageY(:); ...
        upperLeftX(:), upperLeftY(:); ...
        lowerLeakageX(:), lowerLeakageY(:); ...
        leftChannelX(:), leftChannelY(:); ...
        transitionX1(:), transitionY1(:); ...
        transitionX2(:), transitionY2(:); ...
        transitionX3(:), transitionY3(:); ...
        transitionX4(:), transitionY4(:)], ...
        'rows');

    G = triangleGrid(gridPoints);
    G = computeGeometry(G);

    cellCentroids = G.cells.centroids;

    cellsToRemove = ...
        (cellCentroids(:, 1) < 100 - 0.3 | ...
        cellCentroids(:, 1) > 100 + 0.3) & ...
        (cellCentroids(:, 2) < 130 & ...
        cellCentroids(:, 2) > 30);

    G = removeCells(G, cellsToRemove);
    G = makeLayeredGrid(pebi(G), 1);
    G = computeGeometry(G);

    cellCentroids = G.cells.centroids;

    cellsToRemove = ...
        (cellCentroids(:, 1) < 99.9 | ...
        cellCentroids(:, 1) > 100.1) & ...
        (cellCentroids(:, 2) < 130 & ...
        cellCentroids(:, 2) > 30);

    G = removeCells(G, cellsToRemove);
    G = computeGeometry(G);

    cellCentroids = G.cells.centroids;

    cellsToRemove = ...
        cellCentroids(:, 1) < -0.5 - eps & ...
        cellCentroids(:, 2) > 130;

    G = removeCells(G, cellsToRemove);
    G = computeGeometry(G);

    cellCentroids = G.cells.centroids;

    cellsToRemove = ...
        cellCentroids(:, 1) < 99.8 & ...
        cellCentroids(:, 2) < 30;

    G = removeCells(G, cellsToRemove);
    G = computeGeometry(G);
end

cellCentroids = G.cells.centroids;
cellTemplate = ones(G.cells.num, 1);

% Preserve the original variable names used by the plotting sections.
c = cellCentroids;
C = cellTemplate;

%% Rock properties

K0 = 2e-14 .* cellTemplate; % Aquifer permeability, m^2

activeCellMap = G.cells.indexMap;

leakageCells = activeCellMap( ...
    cellCentroids(:, 1) > 99.9 & ...
    cellCentroids(:, 1) < 100.1 & ...
    cellCentroids(:, 2) < 130 & ...
    cellCentroids(:, 2) > 30);

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
     15 * hour, dt, 6e-3, 0.01, 0,    0; ...
     11 * hour, dt, 6e-3, 0,    0,    0; ...
     74 * hour, dt, 0,    0,    0,    0; ...
     30 * hour, dt, 6e-3, 0,    0.04, 0; ...
      5 * hour, dt, 6e-3, 0,    0,    0; ...
     25 * hour, dt, 0,    0,    0,    0; ...
     40 * hour, dt, 6e-3, 0,    0,  300; ...
     10 * hour, dt, 6e-3, 0,    0,    0; ...
    390 * hour, dt, 0,    0,    0,    0; ...
     30 * hour, dt, 6e-3, 0,    0.04, 0; ...
     20 * hour, dt, 6e-3, 0,    0,    0; ...
     20 * hour, dt, 0,    0,    0,    0; ...
     20 * hour, dt, 6e-3, 0,    0,  300; ...
     20 * hour, dt, 6e-3, 0,    0,    0; ...
     90 * hour, dt, 0,    0,    0,    0; ...
     20 * hour, dt, 6e-3, 0,    0,  300; ...
     20 * hour, dt, 6e-3, 0,    0,    0; ...
    110 * hour, dt, 0,    0,    0,    0];

if exist('AD_MICP_TEST', 'var')
    dt = 30 * minute;
    M(:, 1) = dt;
    M(:, 2) = dt;
end

N = size(M, 1);

%% MICP injection well

r = 0.15;
Whu = 1 / 10;
Whb = 1 - Whu;

allCellIndices = (1 : G.cells.num)';

cellsWu = allCellIndices( ...
    cellCentroids(:, 1) < min(cellCentroids(:, 1)) + 0.1 & ...
    cellCentroids(:, 2) > 130 & ...
    cellCentroids(:, 2) < 133);

W = addWell([], G, rock, cellsWu, 'Type', 'rate', ...
    'Comp_i', [1, 0], 'Val', Whu .* M(1, 3), ...
    'Radius', r, 'dir', 'y');

cellsWb = allCellIndices( ...
    cellCentroids(:, 1) < min(cellCentroids(:, 1)) + 0.1 & ...
    cellCentroids(:, 2) > 133);

W = addWell(W, G, rock, cellsWb, 'Type', 'rate', ...
    'Comp_i', [1, 0], 'Val', Whb .* M(1, 3), ...
    'Radius', r, 'dir', 'y');

for wellIndex = 1 : numel(W)
    W(wellIndex).m = 0;
    W(wellIndex).o = 0;
    W(wellIndex).u = 0;
end

W(1).m = M(1, 4);
W(1).o = M(1, 5);
W(1).u = M(1, 6);

G.injectionwellonboundary = 1;
G.cellsinjectionwell = [cellsWu; cellsWb];

%% Gravity

gravity reset
gravity on
gravity y

gravityMagnitude = norm(gravity);

%% Boundary conditions

boundaryFaceIndices = boundaryFaces(G);

openBoundaryFaces = boundaryFaceIndices( ...
    abs(G.faces.normals(boundaryFaceIndices, 1)) > eps & ...
    (G.faces.centroids(boundaryFaceIndices, 1) < -L + 2 | ...
    G.faces.centroids(boundaryFaceIndices, 1) > L - 2));

boundaryPressure = ...
    G.faces.centroids(openBoundaryFaces, 2) .* ...
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
    cellCentroids(:, 2) .* fluid.rhoWS .* gravityMagnitude;

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
    micpPlotFunction = getPlotAfterStepMICP(state0, model, 0, 270);

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
    ntco2 = 10;
end

co2Timesteps = repmat(co2TimeStep, ntco2, 1);

%% CO2 injection well

QCO2 = (1600 / day) / L;

co2WellCells = allCellIndices( ...
    cellCentroids(:, 1) < min(cellCentroids(:, 1)) + 0.1 & ...
    cellCentroids(:, 2) > 130);

co2Well = addWell([], G, rock, co2WellCells, 'Type', 'rate', ...
    'Comp_i', [0, 1], 'Val', QCO2, 'Radius', r, 'dir', 'y');

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

%% Simulate CO2 migration

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    [~, statese] = simulateScheduleAD(co2InitialState, modela, co2Schedule);
    [~, statesf] = simulateScheduleAD(co2InitialState, modelb, co2Schedule);
    [~, statesg] = simulateScheduleAD(co2InitialState, modelc, co2Schedule);
    [~, statesh] = simulateScheduleAD(co2InitialState, modeld, co2Schedule);
else
    co2PlotFunctionA = getPlotAfterStepCO2( ...
        co2InitialState, modela, 0, 270);
    co2PlotFunctionB = getPlotAfterStepCO2( ...
        co2InitialState, modelb, 0, 270);
    co2PlotFunctionC = getPlotAfterStepCO2( ...
        co2InitialState, modelc, 0, 270);
    co2PlotFunctionD = getPlotAfterStepCO2( ...
        co2InitialState, modeld, 0, 270);

    [~, statese] = simulateScheduleAD(co2InitialState, modela, ...
                            co2Schedule, 'afterStepFn', co2PlotFunctionA);
    [~, statesf] = simulateScheduleAD(co2InitialState, modelb, ...
                            co2Schedule, 'afterStepFn', co2PlotFunctionB);
    [~, statesg] = simulateScheduleAD(co2InitialState, modelc, ...
                            co2Schedule, 'afterStepFn', co2PlotFunctionC);
    [~, statesh] = simulateScheduleAD(co2InitialState, modeld, ...
                            co2Schedule, 'afterStepFn', co2PlotFunctionD);
end

%% Compute CO2 leakage rates

allFaceIndices = (1 : G.faces.num)';

leakageFaces = allFaceIndices( ...
    G.faces.centroids(:, 2) < 80.6 & ...
    G.faces.centroids(:, 2) > 80.3 & ...
    abs(G.faces.normals(:, 2)) > 0.1);

lr0 = zeros(ntco2, 1);
lr1 = zeros(ntco2, 1);
lr2 = zeros(ntco2, 1);
lr3 = zeros(ntco2, 1);

for stepIndex = 1 : ntco2
    lr0(stepIndex) = abs(statese{stepIndex}.flux(leakageFaces(1), 2));
    lr1(stepIndex) = abs(statesf{stepIndex}.flux(leakageFaces(1), 2));
    lr2(stepIndex) = abs(statesg{stepIndex}.flux(leakageFaces(1), 2));
    lr3(stepIndex) = abs(statesh{stepIndex}.flux(leakageFaces(1), 2));
end

statese = statese{end};
statesf = statesf{end};
statesg = statesg{end};
statesh = statesh{end};

% Preserve the original timestep variable used in the Figure 11b plotting
% section.
dt = co2TimeStep;

%% Export results in GNU Octave

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    outputDirectory = fullfile(pwd, 'vtk_micp_2Dfls');

    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    micpOutput = fullfile(outputDirectory, 'states_2Dfls');
    co2Output = fullfile(outputDirectory, 'states_2Dfls_CO2');

    mrsttovtk(G, states, micpOutput, '%f');
    mrsttovtk(G, statesh, co2Output, '%f');

    fprintf(['VTK results written to:\n' ...
             '  %s.pvd\n' ...
             '  %s.pvd\n'], ...
             micpOutput, co2Output);

    return
end

% Figure 10 paper (MATLAB)
porosityf = porosity - statesb.c - statesb.b;
porosityg = porosity - statesc.c - statesc.b;
porosityh = porosity - statesd.c - statesd.b;
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
colormap (n1, cc);
caxis([2e-14 1e-12]);
cb = colorbar;
title(cb, '$m^2$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.265 0.67 0.005 0.15], 'Ticks', [2e-14 1e-12], ...
                                                           'FontSize', fS);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]'; '(a)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
plotCellData(G, K0);
title('Initial permeability', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', 0 : 100 : L, 'YTick', 0 : 20 : H, ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 270);
line([200 100], [60 30], [0 0],'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([200 100], [105 130], [0 0], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
ax1 = axes('position', [0.195 0.71 0.045 0.09], 'YAxisLocation', 'right');
box on
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s = plotCellData(G, K0);
s.EdgeColor = 'none';
colormap (ax1, cc);
set(gca, 'FontSize', 6, 'XTick', [99.85 100.15], 'YTick', [30 80 130], ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 270)
n2 = subplot(2, 4, 2);
view(0, 0);
colormap (n2, ccc);
caxis([0 100]);
cb = colorbar;
title(cb, '\%', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.475 0.67 0.005 0.15], 'YTick', [0 100]);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]' ; '(b)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
s = plotCellData(G, 100 * (K0 - fluid.K(porosity - statesb.c - ...
                                                        statesb.b)) ./ K0);
s.EdgeColor = 'none';
title('Permeability (phase I MICP)', 'FontSize', fS, ...
                              'FontName', 'Arial', 'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', 0 : 100 : L, 'YTick', 0 : 20 : H, ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 270);
line([200 100], [60 30], [0 0], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([200 100], [105 130], [0 0], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
ax2 = axes('position', [0.4 0.71 0.045 0.09], 'YAxisLocation', 'right');
box on
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s = plotCellData(G, 100 * (K0 - fluid.K(porosity - statesb.c - ...
                                                        statesb.b)) ./ K0);
s.EdgeColor = 'none';
set(gca, 'FontSize', 6, 'XTick', [99.85 100.15], 'YTick', [30 80 130], ...
                                     'color', 'none', 'FontName', 'Arial');
colormap (ax2, ccc);
caxis([0 100]);
view(0, 270)
n3 = subplot(2, 4, 3);
view(0, 0);
colormap (n3, ccc);
caxis([0 100]);
cb = colorbar;
title(cb, '\%', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.68 0.67 0.005 0.15], 'YTick', [0 100]);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]' ; '(c)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
s = plotCellData(G, 100 * (K0 - fluid.K(porosity - statesc.c - ...
                                                        statesc.b)) ./ K0);
s.EdgeColor = 'none';
title('Permeability (phase II MICP)', 'FontSize', fS, ...
                              'FontName', 'Arial', 'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', 0 : 100 : L, 'YTick', 0 : 20 : H, ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 270);
line([200 100], [60 30], [0 0], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([200 100], [105 130], [0 0], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
ax3 = axes('position', [0.605 0.71 0.045 0.09], 'YAxisLocation', 'right');
box on
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s = plotCellData(G, 100 * (K0 - fluid.K(porosity - statesc.c - ...
                                                        statesc.b)) ./ K0);
s.EdgeColor = 'none';
set(gca, 'FontSize', 6, 'XTick', [99.85 100.15], 'YTick', [30 80 130], ...
                                     'color', 'none', 'FontName', 'Arial');
colormap (ax3, ccc);
caxis([0 100]);
view(0, 270)
n4 = subplot(2, 4, 4);
view(0, 0);
colormap (n4, ccc);
caxis([0 100]);
cb = colorbar;
title(cb, '\%', 'FontSize', fS, 'Interpreter', 'latex', 'FontName', ...
                                                                  'Arial');
set(cb, 'position', [0.89 0.67 0.005 0.15], 'YTick', [0 100]);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]' ; '(d)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
s = plotCellData(G, 100 * (K0 - fluid.K(porosity - statesd.c - ...
                                                        statesd.b)) ./ K0);
s.EdgeColor = 'none';
title('Permeability (phase III MICP)', 'FontSize', fS, ...
                              'FontName', 'Arial', 'Interpreter', 'latex');
set(gca, 'FontSize', fS,'XTick', 0 : 100 : L, 'YTick', 0 : 20 : H, ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 270);
line([200 100], [60 30], [0 0], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([200 100], [105 130], [0 0], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
ax4 = axes('position', [0.815 0.71 0.045 0.09], 'YAxisLocation', 'right');
box on
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s = plotCellData(G, 100 * (K0 - fluid.K(porosity - statesd.c - ...
                                                        statesd.b)) ./ K0);
s.EdgeColor = 'none';
set(gca, 'FontSize', 6, 'XTick', [99.85 100.15], 'YTick', [30 80 130], ...
                                     'color', 'none', 'FontName', 'Arial');
colormap (ax4, ccc);
caxis([0 100]);
view(0, 270)
n5 = subplot(2, 4, 5);
view(0, 0);
colormap (n5, c);
caxis([0 75]);
cb = colorbar;
title(cb, 'kg/m$^3$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.265 0.2 0.005 0.15], 'YTick', [0 25 50 75]);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]' ; '(e)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
s = plotCellData(G, fluid.rhoOS * porosity * statese.s(:, 2));
s.EdgeColor = 'none';
title('CO$_2$ (100 days)', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', 0 : 100 : L, 'YTick', 0 : 20 : H, ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 270);
line([200 100], [60 30], [0 0], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([200 100], [105 130], [0 0], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
ax5 = axes('position', [0.195 0.237 0.045 0.09], 'YAxisLocation', 'right');
box on
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s = plotCellData(G, fluid.rhoOS * porosity * statese.s(:, 2));
s.EdgeColor = 'none';
colormap (ax5, c);
caxis([0 75]);
set(gca,'FontSize', 6, 'XTick', [99.85 100.15], 'YTick', [30 80 130], ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 270)
n6 = subplot(2, 4, 6);
view(0, 0);
colormap (n6, c);
caxis([0 75]);
cb = colorbar;
title(cb, 'kg/m$^3$', 'FontSize', fS, 'Interpreter', ...
                                             'latex', 'FontName', 'Arial');
set(cb, 'position', [0.475 0.2 0.005 0.15], 'YTick', [0 25 50 75]);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]' ; '(f)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
s = plotCellData(G, fluid.rhoOS * porosityf .* statesf.s(:, 2));
s.EdgeColor = 'none';
title('CO$_2$ (phase I MICP)', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', 0 : 100 : L, 'YTick', 0 : 20 : H, ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 270);
line([200 100], [60 30], [0 0], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([200 100], [105 130], [0 0],'Color','black','LineStyle','--', ...
                                                          'LineWidth', lW);
ax6 = axes('position', [0.4 0.237 0.045 0.09], 'YAxisLocation', 'right');
box on
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s = plotCellData(G, fluid.rhoOS * porosityf .* statesf.s(:, 2));
s.EdgeColor = 'none';
colormap (ax6, c);
caxis([0 75]);
set(gca, 'FontSize', 6, 'XTick', [99.85 100.15], 'YTick', [30 80 130], ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 270)
n7 = subplot(2, 4, 7);
view(0, 0);
colormap(n7, c);
caxis([0 75]);
cb = colorbar;
title(cb, 'kg/m$^3$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.68 0.2 0.005 0.15], 'YTick', [0 25 50 75]);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]' ; '(g)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
s = plotCellData(G, fluid.rhoOS * porosityg .* statesg.s(:, 2));
s.EdgeColor = 'none';
title('CO$_2$ (phase II MICP)', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', 0 : 100 : L, 'YTick', 0 : 20 : H, ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 270);
line([200 100], [60 30], [0 0], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([200 100], [105 130], [0 0], 'Color', 'black','LineStyle','--', ...
                                                          'LineWidth', lW);
ax7 = axes('position', [0.605 0.237 0.045 0.09], 'YAxisLocation', 'right');
box on
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s = plotCellData(G,fluid.rhoOS * porosityg .* statesg.s(:, 2));
s.EdgeColor = 'none';
colormap (ax7, c);
caxis([0 75]);
set(gca, 'FontSize', 6, 'XTick', [99.85 100.15], 'YTick', [30 80 130], ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 270)
n8 = subplot(2, 4, 8);
view(0, 0);
colormap (n8, c);
caxis([0 75]);
cb = colorbar;
title(cb, 'kg/m$^3$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.89 0.2 0.005 0.15], 'YTick', [0 25 50 75]);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]' ; '(h)'}, 'FontSize', fS, 'FontName', 'Arial');
ylabel('z [m]', 'FontSize', fS, 'FontName', 'Arial');
s = plotCellData(G, fluid.rhoOS * porosityh .* statesh.s(:, 2));
s.EdgeColor = 'none';
title('CO$_2$ (phase III MICP)', 'FontSize', fS, 'FontName', ...
                                          'Arial', 'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'XTick', 0 : 100 : L, 'YTick', 0 : 20 : H, ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 270);
line([200 100], [60 30], [0 0], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
line([200 100], [105 130], [0 0], 'Color', 'black', 'LineStyle', '--', ...
                                                          'LineWidth', lW);
ax8 = axes('position', [0.815 0.237 0.045 0.09], 'YAxisLocation', 'right');
box on
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s = plotCellData(G, fluid.rhoOS * porosityh .* statesh.s(:, 2));
s.EdgeColor = 'none';
colormap (ax8, c);
caxis([0 75]);
set(gca, 'FontSize', 6, 'XTick', [99.85 100.15], 'YTick', [30 80 130], ...
                                     'color', 'none', 'FontName', 'Arial');
view(0, 270)
%print -depsc2 Fig10.eps

% Figures 11a and 11b paper
clear m o u b c K vc v cell_leak
cells = 1 : G.cells.num;
nt_micp = cumsum(schedule.step.val) / hour;

cell_leak = cells(G.cells.centroids(:, 2) < 130 & ...
                                             G.cells.centroids(:, 2) > 30);
for i = 1 : nt
  c(i) = mean(states{i}.c(cell_leak));
  b(i) = mean(states{i}.b(cell_leak));
  m(i) = mean(states{i}.m(cell_leak));
  u(i) = mean(states{i}.u(cell_leak));
  o(i) = mean(states{i}.o(cell_leak));
  currentPorosity = max(porosity - states{i}.c - states{i}.b, ...
                      model.minimumPorosity);
  Ki = fluid.K(currentPorosity);
  K(i) = mean(Ki(cell_leak) ./ K0(cell_leak));
  vc = faceFlux2cellVelocity(G, states{i}.flux(:));
  v(i) = mean(sqrt(sum(vc(cell_leak, :) .^ 2, 2)));
end
fS = 11;
lW = 3;
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
text(650, 1.02, 'Phase II', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
text(825, 1.02, 'Phase III', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
xlim([0 nt_micp(end)]);
xlabel({'Time [h]' ; '(a)'}, 'FontSize', fS, 'Interpreter', 'latex');
ylabel('[$-$]', 'FontSize', fS, 'Interpreter', 'latex');
h = legend('$v_w/0.0070\textrm{ m/s}$', '$c_m/0.0020\textrm{ kg/m}^3$', ...
          '$c_o/0.0084\textrm{ kg/m}^3$', '$c_u/136 \textrm{ kg/m}^3$', ...
                                  '$\phi_b/0.0003$', '$\phi_c/0.0340$', ...
                                           '$K/10^{-12}\textrm{ m}^2$', ...
                                                'Interpreter', 'latex', ...
                                                           'FontSize', fS);
rect = [0.37, 0.35, 0.2, 0.25];
set(h, 'Position', rect);
set(gca, 'FontSize', fS, 'FontName', 'Arial', 'XTick', 0 : 100 : 1000, ...
                                             'YGrid', 'on', 'XGrid', 'on');
%print -depsc2 Fig11a.eps

fS = 11;
lW = 9;
figure('Units', 'inches', 'PaperPosition', [0 0 6.83 6.83]);
hold on
plot((1 : ntco2) * dt / day, 100 * lr0 / QCO2, ...
                  'color', [1 0.2 0.2], 'LineWidth', lW, 'LineStyle', '-');
plot((1 : ntco2) * dt / day, 100 * lr1 / QCO2, ...
                    'color', [1 0.5 0], 'LineWidth', lW, 'LineStyle', '-');
plot((1 : ntco2) * dt / day, 100 * lr2 / QCO2, ...
             'color', [0.61 0.61 0.61], 'LineWidth', lW, 'LineStyle', '-');
plot((1 : ntco2) * dt / day, 100 * lr3 / QCO2, ...
                        'color',[0 0 0], 'LineWidth', lW, 'LineStyle','-');
hold off
xlim([0 100]);
ylim([0 60]);
xlabel({'Time [d]' ; '(b)'}, 'FontSize', fS, 'Interpreter', 'latex');
ylabel('CO$_2$ leakage rate/injection rate [\%]', 'FontSize', fS, ...
                                                   'Interpreter', 'latex');
grid on
legend('Without MICP', 'Phase I MICP', 'Phase II MICP', ...
                                     'Phase III MICP', 'Location', 'best');
set(gca, 'FontSize', fS, 'FontName', 'Arial', 'XTick', 0 : 20 : 100, ...
                                                     'YTick', 0 : 10 : 60);
%print -depsc2 Fig11b.eps 
