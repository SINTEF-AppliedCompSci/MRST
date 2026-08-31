%% One-dimensional horizontal MICP treatment
%
% Set up and solve the one-dimensional horizontal-flow system (1Dfhs).
%
% In MATLAB, this example produces Figure 5 from:
%
% Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E., 2021. Practical
% approaches to study microbially induced calcite precipitation at the
% field scale. International Journal of Greenhouse Gas Control 106, 103256.
%
% https://doi.org/10.1016/j.ijggc.2021.103256
%
% In GNU Octave, the example writes the results to the
% `vtk_micp_1Dfhs` directory for visualization in ParaView.
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

L = 75;   % Aquifer length, m
l = 25;   % Region where MICP processes are most relevant, m
dw = 0.5; % Size of the well cell, m
dl = 0.05;% Size of the cells within the treatment region, m

X = [0, dw, dw + dl : dl : l, L .* exp(-1.075 : 0.025 : 0)];
X(end) = L;

G = tensorGrid(X, [0 1], [0 1]);
G = computeGeometry(G);

cellTemplate = ones(G.cells.num, 1);

%% Rock properties

K0 = 1e-12 .* cellTemplate; % Initial permeability, m^2
porosity = 0.2;             % Initial porosity, [-]

rock = makeRock(G, K0, porosity);

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
    (porosity - fluid.crit), 0);

fluid.K = @(currentPorosity) max( ...
    (K0 .* normalizedPorosity(currentPorosity) .^ fluid.eta + ...
    fluid.kmin) .* K0 ./ (K0 + fluid.kmin), fluid.kmin);

%% Injection strategy
%
% Each row of `M` contains:
%
%   1. Phase duration, s
%   2. Simulation timestep, s
%   3. Injection rate, m^3/s
%   4. Microorganism concentration, kg/m^3
%   5. Oxygen concentration, kg/m^3
%   6. Urea concentration, kg/m^3

dt_on = 20 * minute;
dt_off = 5 * hour;

M = [ ...
     20 * hour, dt_on,  2.4e-5, 0.01, 0,    0; ...
     20 * hour, dt_on,  2.4e-5, 0,    0,    0; ...
    100 * hour, dt_off, 0,      0,    0,    0; ...
     20 * hour, dt_on,  2.4e-5, 0,    0.04, 0; ...
     20 * hour, dt_on,  2.4e-5, 0,    0,    0; ...
     50 * hour, dt_off, 0,      0,    0,    0; ...
     20 * hour, dt_on,  2.4e-5, 0,    0,  300; ...
     20 * hour, dt_on,  2.4e-5, 0,    0,    0; ...
    230 * hour, dt_off, 0,      0,    0,    0];

if exist('AD_MICP_TEST', 'var')
    dt_on = 30 * minute;
    M = [5 * hour, dt_on, 2.4e-5, 0.01, 0.04, 300];
end

N = size(M, 1);

%% Injection well

r = 0.15;

W = addWell([], G, rock, 1, 'Type', 'rate', 'Comp_i', [1, 0], ...
                                             'Val', M(1, 3), 'Radius', r);

W.m = M(1, 4);
W.o = M(1, 5);
W.u = M(1, 6);

% The injection well is located on the left boundary. These fields are used
% to correct the reconstructed cell velocity when computing dispersion and
% detachment.
G.injectionwellonboundary = 1;
G.cellsinjectionwell = 1;

%% Outflow boundary condition

boundaryFaceIndices = boundaryFaces(G);

outflowFaces = boundaryFaceIndices( ...
    abs(G.faces.normals(boundaryFaceIndices, 1)) > eps & ...
    G.faces.centroids(boundaryFaceIndices, 1) > X(end - 1));

bc = addBC([], outflowFaces, 'pressure', atm, 'sat', [0 0]);

numberOfBoundaryFaces = size(bc.sat, 1);

bc.m = zeros(numberOfBoundaryFaces, 1);
bc.o = zeros(numberOfBoundaryFaces, 1);
bc.u = zeros(numberOfBoundaryFaces, 1);
bc.b = zeros(numberOfBoundaryFaces, 1);
bc.c = zeros(numberOfBoundaryFaces, 1);

%% Construct simulation schedule

stepsPerPhase = round(M(:, 1) ./ M(:, 2));
totalNumberOfSteps = sum(stepsPerPhase);

timesteps = zeros(totalNumberOfSteps, 1);
controlIndices = zeros(totalNumberOfSteps, 1);

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
    schedule.control(phaseIndex).W.val = M(phaseIndex, 3);
    schedule.control(phaseIndex).W.m = M(phaseIndex, 4);
    schedule.control(phaseIndex).W.o = M(phaseIndex, 5);
    schedule.control(phaseIndex).W.u = M(phaseIndex, 6);
end

schedule.step.control = controlIndices;
schedule.step.val = timesteps;

phaseEndSteps = cumsum(stepsPerPhase);

% Maximum injected oxygen and urea concentrations
fluid.Comax = max(M(:, 5));
fluid.Cumax = max(M(:, 6));

%% Create model and solver

model = MICPModel(G, rock, fluid);
model.toleranceMB = 1e-14;
model.nonlinearTolerance = 1e-14;

solver = getNonLinearSolver(model);
solver.LinearSolver.tolerance = 1e-14;

%% Initial state

state0 = initState(G, W, atm, [1, 0]);

state0.m = zeros(G.cells.num, 1);
state0.o = zeros(G.cells.num, 1);
state0.u = zeros(G.cells.num, 1);
state0.b = zeros(G.cells.num, 1);
state0.c = zeros(G.cells.num, 1);

%% Run the simulation

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    [~, states] = simulateScheduleAD(state0, model, schedule, ...
                                     'NonLinearSolver', solver);
else
    afterStepFunction = getPlotAfterStepMICP(state0, model, 0, 90);

    [~, states] = simulateScheduleAD(state0, model, schedule, ...
                             'NonLinearSolver', solver, ...
                             'afterStepFn', afterStepFunction);
end

%% Export results in GNU Octave

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    outputDirectory = fullfile(pwd, 'vtk_micp_1Dfhs');

    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    outputName = fullfile(outputDirectory, 'states');
    mrsttovtk(G, states, outputName, '%f');

    fprintf('VTK results written to:\n  %s.pvd\n', outputName);
    return
end

%% Figure 5

lW = 2;
fS = 9;
pL = [10 15];
lS = {'-', ':', ':', '-', '-.', '-.', '-', '--', '--'};
pltCls = {[0 0.8 0], [0 0.74 1], [0 0 0], [1 0.5 0.9], [0 0.74 1], ...
                                  [0 0 0], [1 0.9 0], [0 0.74 1], [0 0 0]};

figure;
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6.83 6]);

n1 = subplot(3, 3, 4);
hold on
for i = 1 : N
    plot(G.cells.centroids(:, 1), states{phaseEndSteps(i)}.m, ...
        'color', pltCls{i}, 'LineWidth', lW, 'LineStyle', lS{i});
end
line([pL(1) pL(1)], [0 W(1).m], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(2) pL(2)], [0 W(1).m], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [W(1).m W(1).m], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [0 0], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
hold off
xlim([0 L]);
ylim([0 W(1).m]);
xlabel({'x [m]' ; '(a)'}, 'FontSize', fS, 'Interpreter', 'latex');
ylabel('$c_m$ [kg/m$^3$]', 'FontSize', fS, 'Interpreter', 'latex');
grid on
title('Microorganisms', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'FontName', 'Arial', ...
                       'XTick', (0 : 10 : L), 'YTick', (0 : 0.002 : 0.01));

n2 = subplot(3, 3, 5);
hold on
for i = 1 : N
    plot(G.cells.centroids(:, 1), states{phaseEndSteps(i)}.o, ...
        'color', pltCls{i}, 'LineWidth', lW, 'LineStyle', lS{i});
end
line([pL(1) pL(1)], [0 fluid.Comax], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(2) pL(2)], [0 fluid.Comax], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [fluid.Comax fluid.Comax], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [0 0], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
hold off
xlim([0 L]);
ylim([0 0.04]);
xlabel({'x [m]' ; '(b)'}, 'FontSize', fS, 'Interpreter', 'latex');
ylabel('$c_o$ [kg/m$^3$]', 'FontSize', fS, 'Interpreter', 'latex');
grid on
title('Oxygen', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'FontName', 'Arial', ...
                        'XTick', (0 : 10 : L), 'YTick', (0 : 0.01 : 0.04));

n3 = subplot(3, 3, 6);
hold on
for i = 1 : N
    plot(G.cells.centroids(:, 1), states{phaseEndSteps(i)}.u, ...
        'color', pltCls{i}, 'LineWidth', lW, 'LineStyle', lS{i});
end
line([pL(1) pL(1)], [0 fluid.Cumax], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(2) pL(2)], [0 fluid.Cumax], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [fluid.Cumax fluid.Cumax], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [0 0], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
hold off
xlim([0 L]);
ylim([0 300]);
xlabel({'x [m]' ; '(c)'}, 'FontSize', fS, 'Interpreter', 'latex');
ylabel('$c_u$ [kg/m$^3$]', 'FontSize', fS, 'Interpreter', 'latex');
grid on
title('Urea', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'FontName', 'Arial', ...
                           'XTick', (0 : 10 : L), 'YTick', (0 : 60 : 300));

n4 = subplot(3, 3, 7);
hold on
for i = 1 : N
    plot(G.cells.centroids(:, 1), states{phaseEndSteps(i)}.b, ...
        'color', pltCls{i}, 'LineWidth', lW, 'LineStyle', lS{i});
end
line([pL(1) pL(1)], [0 0.0003], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(2) pL(2)], [0 0.0003], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [0.0003 0.0003], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [0 0], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
hold off
xlim([0 L]);
ylim([0 0.0003]);
xlabel({'x [m]' ; '(d)'}, 'FontSize', fS, 'Interpreter', 'latex');
ylabel('$\phi_b$ [$-$]', 'FontSize', fS, 'Interpreter', 'latex');
grid on
title('Biofilm', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'FontName', 'Arial', ...
                    'XTick', (0 : 10 : L), 'YTick', (0 : 0.0001 : 0.0003));

n5 = subplot(3, 3, 8);
hold on
for i = 1 : N
    plot(G.cells.centroids(:, 1), states{phaseEndSteps(i)}.c, ...
        'color', pltCls{i}, 'LineWidth', lW, 'LineStyle', lS{i});
end
line([pL(1) pL(1)], [0 0.04], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(2) pL(2)], [0 0.04], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [0.04 0.04], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [0 0], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
hold off
xlim([0 L]);
ylim([0 0.04]);
xlabel({'x [m]' ; '(e)'}, 'FontSize', fS, 'Interpreter', 'latex');
ylabel('$\phi_c$ [$-$]', 'FontSize', fS, 'Interpreter', 'latex');
grid on
title('Calcite', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'FontName', 'Arial', ...
                         'XTick', (0 : 10 : L), 'YTick', (0 : 0.01 : 0.05));

n6 = subplot(3, 3, 9);
hold on
for i = 1 : N
    currentPorosity = max(porosity - states{phaseEndSteps(i)}.b - ...
                              states{phaseEndSteps(i)}.c, model.minimumPorosity);
    plot(G.cells.centroids(:, 1), 100 .* (1 - fluid.K(currentPorosity) ./ K0), ...
        'color', pltCls{i}, 'LineWidth', lW, 'LineStyle', lS{i});
end
line([pL(1) pL(1)], [0 100], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(2) pL(2)], [0 100], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [100 100], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [0 0], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
hold off
xlim([0 L]);
ylim([0 100]);
xlabel({'x [m]' ; '(f)'}, 'FontSize', fS, 'Interpreter', 'latex');
ylabel('$|\Delta K/K_0|$ [$\%$]', 'FontSize', fS, 'Interpreter', 'latex');
grid on
title('Permeability', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');

cb = legend('$t^I_1=\;\;20\;$h', '$t^I_2=\;\;40\;$h', ...
               '$t^I_3=140\;$h$\qquad \qquad \qquad$', '$t^I_4=160\;$h', ...
               '$t^I_5=180\;$h', '$t^I_6=230\;$h$\qquad \qquad \qquad$', ...
                    '$t^I_7=250\;$h', '$t^I_8=270\;$h', '$t^I_9=500\;$h', ...
                   'Location', 'best', 'Interpreter', 'latex', 'FontSize', fS);

set(cb, 'position', [0.5 0.67 0.01 0.15]);

if isprop(cb, 'NumColumns')
    cb.NumColumns = 3;
end

set(gca, 'FontSize', fS, 'FontName', 'Arial', ...
                             'XTick', (0 : 10 : L), 'YTick', (0 : 20 : 100));

% print -depsc2 Fig5.eps
