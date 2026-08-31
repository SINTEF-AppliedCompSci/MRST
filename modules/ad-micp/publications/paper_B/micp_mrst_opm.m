%% MRST and OPM comparison for the MICP model
%
% This example sets up and solves a one-dimensional MICP treatment problem
% in MRST and compares the results with corresponding OPM simulation data.
%
% In MATLAB, this file produces Figure 3 in:
%
% Landa-Marbán, D., Kumar, K., Tveit, S., Gasda, S.E. Numerical studies of
% CO2 leakage remediation by MICP-based plugging technology. In: Røkke,
% N.A. and Knuutila, H.K. (Eds.), Short Papers from the 11th International
% Trondheim CCS Conference, ISBN 978-82-536-1714-5, pp. 284-290.
%
% In GNU Octave, this file writes the simulation results to the
% `vtk_micp_mrst_opm` folder for visualization in ParaView.
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

% Required modules
mrstModule add ad-blackoil ad-core ad-micp

%% Grid
domainLength = 100;  % Aquifer length, m
G = tensorGrid(0 : domainLength, [0 1], [0 1]);
G = computeGeometry(G);
cellTemplate = ones(G.cells.num, 1);

%% Rock properties
initialPermeability = 1e-14 .* cellTemplate;  % Permeability, m^2
initialPorosity = 0.15;                       % Porosity, [-]
rock = makeRock(G, initialPermeability, initialPorosity);

%% Fluid and MICP model properties
fluid.muw = 2.535e-4;        % Water viscosity, Pa s
fluid.bW = @(pressure) ...
    0 .* pressure + 1;       % Water formation volume factor, [-]
fluid.bO = fluid.bW;         % Compatibility field, [-]
fluid.rhoWS = 1045;          % Water density, kg/m^3
fluid.rhoOS = 479;           % CO2 density, kg/m^3

fluid.rho_b = 35;            % Biofilm density, kg/m^3
fluid.rho_c = 2710;          % Calcite density, kg/m^3
fluid.k_str = 2.6e-10;       % Detachment coefficient, m/(Pa s)
fluid.diffm = 0;             % Microorganism diffusion, m^2/s
fluid.diffo = 0;             % Oxygen diffusion, m^2/s
fluid.diffu = 0;             % Urea diffusion, m^2/s
fluid.alphaL = 0;            % Longitudinal dispersivity, m
fluid.alphaT = 0;            % Transverse dispersivity, m
fluid.eta = 3;               % Permeability fitting factor, [-]
fluid.k_o = 2e-5;            % Oxygen half-velocity constant, kg/m^3
fluid.k_u = 21.3;            % Urea half-velocity constant, kg/m^3
fluid.mu = 4.17e-5;          % Maximum specific growth rate, 1/s
fluid.mu_u = 0.0161;         % Maximum urease utilization rate, 1/s
fluid.k_a = 8.51e-7;         % Microorganism attachment rate, 1/s
fluid.k_d = 3.18e-7;         % Microorganism death rate, 1/s
fluid.Y = 0.5;               % Growth yield coefficient, [-]
fluid.Yuc = 1.67;            % Calcite-to-urea yield coefficient, [-]
fluid.F = 0.5;               % Oxygen consumption factor, [-]
fluid.crit = 0.1;            % Critical porosity, [-]
fluid.kmin = 1e-20;          % Minimum permeability, m^2

% Porosity-permeability relationship
normalizedPorosity = @(porosity) max( ...
    (porosity - fluid.crit) ./ ...
    (initialPorosity - fluid.crit), 0);
fluid.K = @(porosity) max( ...
    (initialPermeability .* ...
    normalizedPorosity(porosity) .^ fluid.eta + ...
    fluid.kmin) .* initialPermeability ./ ...
    (initialPermeability + fluid.kmin), ...
    fluid.kmin);

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
shutInTimeStep = 5 * hour;
injectionStrategy = [ ...
     20 * hour, activeTimeStep, 2 / day, 0.01, 0,    0; ...
     20 * hour, activeTimeStep, 2 / day, 0,    0,    0; ...
    100 * hour, shutInTimeStep, 0,       0,    0,    0; ...
     20 * hour, activeTimeStep, 2 / day, 0,    0.04, 0; ...
     20 * hour, activeTimeStep, 2 / day, 0,    0,    0; ...
     50 * hour, shutInTimeStep, 0,       0,    0,    0; ...
     20 * hour, activeTimeStep, 2 / day, 0,    0,  300; ...
     20 * hour, activeTimeStep, 2 / day, 0,    0,    0; ...
    230 * hour, shutInTimeStep, 0,       0,    0,    0];

if exist('AD_MICP_TEST', 'var')
    injectionStrategy = [ ...
        10 * hour, activeTimeStep, 2 / day, 0.01, 0.04, 300];
end
numberOfPhases = size(injectionStrategy, 1);

%% Injection and pressure-control wells
wellRadius = 0.15;
W = addWell([], G, rock, 1, ...
    'Type', 'rate', 'Comp_i', [1, 0], ...
    'Val', injectionStrategy(1, 3), 'Radius', wellRadius);
W = addWell(W, G, rock, G.cells.num, ...
    'Type', 'bhp', 'Comp_i', [1, 0], ...
    'Val', atm, 'Radius', wellRadius);
% The injection well lies on the boundary. This information is used to
% correct the reconstructed cell velocity in the dispersion calculation.
G.injectionwellonboundary = 1;
G.cellsinjectionwell = 1;
for wellIndex = 1 : numel(W)
    W(wellIndex).m = 0;
    W(wellIndex).o = 0;
    W(wellIndex).u = 0;
end
W(1).m = injectionStrategy(1, 4);
W(1).o = injectionStrategy(1, 5);
W(1).u = injectionStrategy(1, 6);

%% Construct simulation schedule
stepsPerPhase = round( ...
    injectionStrategy(:, 1) ./ injectionStrategy(:, 2));
totalNumberOfSteps = sum(stepsPerPhase);
timesteps = zeros(totalNumberOfSteps, 1);
controlIndices = zeros(totalNumberOfSteps, 1);
firstStep = 1;

for phaseIndex = 1 : numberOfPhases
    lastStep = ...
        firstStep + stepsPerPhase(phaseIndex) - 1;

    timesteps(firstStep : lastStep) = ...
        injectionStrategy(phaseIndex, 2);

    controlIndices(firstStep : lastStep) = ...
        phaseIndex;

    firstStep = lastStep + 1;
end
schedule = simpleSchedule(timesteps, 'W', W);
baseControl = schedule.control(1);
schedule.control = repmat(baseControl, numberOfPhases, 1);
for phaseIndex = 1 : numberOfPhases
    schedule.control(phaseIndex).W(1).val = ...
        injectionStrategy(phaseIndex, 3);

    schedule.control(phaseIndex).W(1).m = ...
        injectionStrategy(phaseIndex, 4);

    schedule.control(phaseIndex).W(1).o = ...
        injectionStrategy(phaseIndex, 5);

    schedule.control(phaseIndex).W(1).u = ...
        injectionStrategy(phaseIndex, 6);
end
schedule.step.control = controlIndices;
schedule.step.val = timesteps;
phaseEndSteps = cumsum(stepsPerPhase);
% Maximum injected oxygen and urea concentrations
fluid.Comax = max(injectionStrategy(:, 5));
fluid.Cumax = max(injectionStrategy(:, 6));

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
    [~, states] = simulateScheduleAD( ...
        state0, model, schedule, ...
        'NonLinearSolver', solver);
else
    afterStepFunction = ...
        getPlotAfterStepMICP(state0, model, 0, 270);

    [~, states] = simulateScheduleAD( ...
        state0, model, schedule, ...
        'NonLinearSolver', solver, ...
        'afterStepFn', afterStepFunction);
end

%% Export results in GNU Octave
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    outputDirectory = fullfile(pwd, 'vtk_micp_mrst_opm');
    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end
    outputName = fullfile(outputDirectory, 'states');
    mrsttovtk(G, states, outputName, '%f');
    fprintf( ...
        'VTK results written to:\n  %s.pvd\n', ...
        outputName);
    return
end

% Read vtk data to states_OPM from the OPM simulation (MATLAB)
fid = fopen('micp_opm_vtk/MICP_OPM.pvd', 'r');
c1 = fscanf(fid, '%c');
fclose(fid);
newStr = extractBetween(c1, '<DataSet timestep="', '" file');
outNums = str2double(newStr);
states_OPM = deal(cell(size(outNums, 1) - 1, 1));
opmNames = {'bacteria concentration', 'oxygen concentration', ...
                               'urea concentration', 'biofilm', 'calcite'};
mrstNames = {'m', 'o', 'u', 'b', 'c'};                          
for i = 0 : max(size(outNums)) - 1
    filename = sprintf('micp_opm_vtk/MICP_OPM-%05d.vtu', i);
    fileID = fopen(filename, 'r');
    c1 = fscanf(fileID, '%c');
    for j = 1 : 5
        s = ['<DataArray type="Float32" Name' ...
              '="' opmNames{j} '" NumberOfComponents="1" format="ascii">'];
        newStr = extractBetween(c1, s, '</DataArray>');
        states_OPM{i + 1}.(mrstNames{j}) = sscanf(newStr{1}, '%f');
    end
    fclose(fileID);
end

% Figure (MATLAB)
lW = 2;
mS = 5;
fS = 9;
pL = [12.5 17.5];
lS = {'-', ':', ':', '-', '-.', '-.', '-', '--', '--'};
pltCls = {[0 0.8 0], [0 0.74 1], [0 0 0], [1 0.5 0.9], [0 0.74 1], ...
                                  [0 0 0], [1 0.9 0], [0 0.74 1], [0 0 0]};                              
figure
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 16 27]);
n2 = subplot(5, 1, 2);
hold on
plot(0, 0, 'color', [0 0 0], 'LineStyle', '-');
for i = 1 : numberOfPhases 
    plot(G.cells.centroids(1 : 2 : end, 1), ...
            states_OPM{i + 1}.m(1 : 2 : end, 1), 'o', 'MarkerSize', mS, ...
                          'MarkerEdgeColor', [0 0 1], 'LineStyle', 'none');
    plot(G.cells.centroids(:, 1), ...
         states{phaseEndSteps(i)}.m, 'color', pltCls{i}, ...
                                      'LineWidth', lW, 'LineStyle', lS{i});                   
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
xlabel('x [m]', 'FontSize', fS, 'Interpreter', 'latex');        
ylabel('$c_m$ [kg/m$^3$]', 'FontSize', fS, 'Interpreter', 'latex');
ax = gca;
ax.YAxis.Exponent = -2;
grid on
title('(a) Microbes', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
legend('MRST', 'OPM', 'Location', 'northeast');
set(gca, 'FontSize', fS, 'FontName', 'Arial', 'XTick', (0 : 20 : L), ...
                                             'YTick', (0 : 0.002: W(1).m));
n5 = subplot(5, 2, 5);
hold on
plot(0,0, 'color', [0 0 0], 'LineStyle', '-');
for i = 1 : numberOfPhases 
    plot(G.cells.centroids(1 : 2 : end, 1), ...
            states_OPM{i + 1}.o(1 : 2 : end, 1), 'o', 'MarkerSize', mS, ...
                          'MarkerEdgeColor', [0 0 1], 'LineStyle', 'none');
    plot(G.cells.centroids(:, 1), ...
         states{phaseEndSteps(i)}.o, 'color', pltCls{i}, ...
                                      'LineWidth', lW, 'LineStyle', lS{i});                   
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
xlim([0 40]);
ylim([0 fluid.Comax]);
xlabel('x [m]', 'FontSize', fS, 'Interpreter', 'latex');        
ylabel('$c_o$ [kg/m$^3$]', 'FontSize', fS, 'Interpreter', 'latex');
ax = gca;
ax.YAxis.Exponent = -2;
grid on
title('(b) Oxygen', 'FontSize', fS, 'FontName', 'Arial', 'Interpreter', ...
                                                                  'latex');
legend('MRST', 'OPM', 'Location', 'northeast');
set(gca, 'FontSize', fS, 'FontName', 'Arial', 'XTick', (0 : 10 : 40), ...
                                        'YTick', (0 : 0.01 : fluid.Comax));
n6 = subplot(5, 2, 6);
hold on
plot(0, 0, 'color',[0 0 0], 'LineStyle', '-');
for i = 1 : numberOfPhases 
    plot(G.cells.centroids(1 : 2 : end, 1), ...
            states_OPM{i + 1}.b(1 : 2 : end, 1), 'o', 'MarkerSize', mS, ...
                          'MarkerEdgeColor', [0 0 1], 'LineStyle', 'none');
    plot(G.cells.centroids(:, 1), ...
         states{phaseEndSteps(i)}.b, 'color', pltCls{i}, ...
                                      'LineWidth', lW, 'LineStyle', lS{i});                   
end
line([pL(1) pL(1)], [0 0.00015], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(2) pL(2)], [0 0.00015], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [0.00015 0.00015], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [0 0], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
hold off
xlim([0 40]);
ylim([0 0.00015]);
xlabel('x [m]', 'FontSize', fS, 'Interpreter', 'latex');        
ylabel('$\phi_b$ [$-$]', 'FontSize', fS, 'Interpreter', 'latex');
grid on
title('(d) Biofilm', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
legend('MRST', 'OPM', 'Location', 'northeast');
set(gca, 'FontSize', fS, 'FontName', 'Arial', 'XTick', (0 : 10 : 40), ...
                                         'YTick', (0 : 0.00003 : 0.00015));
n7 = subplot(5, 2, 7);
hold on
plot(0, 0, 'color', [0 0 0], 'LineStyle', '-');
for i = 1 : numberOfPhases 
    plot(G.cells.centroids(1 : 2 : end, 1), ...
            states_OPM{i + 1}.u(1 : 2 : end, 1), 'o', 'MarkerSize', mS, ...
                          'MarkerEdgeColor', [0 0 1], 'LineStyle', 'none');
    plot(G.cells.centroids(:, 1), ...
         states{phaseEndSteps(i)}.u, 'color', pltCls{i}, ...
                                      'LineWidth', lW, 'LineStyle', lS{i});                   
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
xlim([0 40]);
ylim([0 fluid.Cumax]);
xlabel('x [m]', 'FontSize', fS, 'Interpreter', 'latex');        
ylabel('$c_u$ [kg/m$^3$]', 'FontSize', fS, 'Interpreter', 'latex');
grid on
title('(c) Urea', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
legend('MRST', 'OPM', 'Location', 'northeast');
set(gca, 'FontSize', fS, 'FontName', 'Arial', 'XTick', (0 : 10 : 40), ...
                                          'YTick', (0 : 60 : fluid.Cumax));                                     
n8 = subplot(5, 2, 8);
hold on
plot(0, 0, 'color', [0 0 0], 'LineStyle', '-');
for i = 1 : numberOfPhases 
    plot(G.cells.centroids(1 : 2 : end, 1), ...
            states_OPM{i + 1}.c(1 : 2 : end, 1), 'o', 'MarkerSize', mS, ...
                          'MarkerEdgeColor', [0 0 1], 'LineStyle', 'none');
    plot(G.cells.centroids(:, 1), ...
         states{phaseEndSteps(i)}.c, 'color', pltCls{i}, ...
                                      'LineWidth', lW, 'LineStyle', lS{i});                   
end
line([pL(1) pL(1)], [0 0.02], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);                                                  
line([pL(2) pL(2)], [0 0.02], 'color','red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [0.02 0.02], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
line([pL(1) pL(2)], [0 0], 'color', 'red', ...
                                        'LineStyle', ':', 'LineWidth', lW);
hold off
xlim([0 40]);
ylim([0 0.02]);
xlabel('x [m]', 'FontSize', fS, 'Interpreter', 'latex');        
ylabel('$\phi_c$ [$-$]', 'FontSize', fS, 'Interpreter', 'latex');
ax = gca;
ax.YAxis.Exponent = -2;
ax.YAxis.TickLabelFormat = '%.1f';
grid on
title('(e) Calcite', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
legend('MRST', 'OPM', 'Location', 'northeast');
set(gca, 'FontSize', fS, 'FontName', 'Arial', 'XTick', (0 : 10 : 40), ...
                                                  'YTick', (0:0.005:0.02));                                                        
n9 = subplot(5, 2, 9);
hold on
plot(0, 0, 'color', [0 0 0], 'LineStyle', '-');
for i = 1 : numberOfPhases 
    plot(G.cells.centroids(1 : 2 : end, 1), 100 / porosity * ( ...
              states_OPM{i + 1}.b(1 : 2 : end, 1) + states_OPM{i + 1}.c ...
                              (1 : 2 : end, 1)), 'o', 'MarkerSize', mS, ...
                               'MarkerEdgeColor', [0 0 1], 'LineStyle', ...
                                                                   'none');
    plot(G.cells.centroids(:, 1), 100 / porosity * ( ...
                            states{phaseEndSteps(i)}.b + ...
                            states{phaseEndSteps(i)}.c), ... 
                                   'color', pltCls{i}, 'LineWidth', lW, ...
                                                       'LineStyle', lS{i});                   
end
line([pL(1) pL(1)], [0 20], 'color', 'red', 'LineStyle', ':', ...
                                                          'LineWidth', lW);                                                  
line([pL(2) pL(2)], [0 20], 'color', 'red', 'LineStyle', ':', ...
                                                          'LineWidth', lW);
line([pL(1) pL(2)], [20 20], 'color', 'red', 'LineStyle', ':', ...
                                                          'LineWidth', lW);
line([pL(1) pL(2)], [0 0], 'color', 'red', 'LineStyle', ':', ...
                                                          'LineWidth', lW);
hold off
xlim([0 40]);
ylim([0 20]);
xlabel('x [m]', 'FontSize', fS, 'Interpreter', 'latex');        
ylabel('$\phi$ [$\%$]', 'FontSize', fS, 'Interpreter', 'latex');
grid on
title('(f) Porosity', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
legend('MRST', 'OPM', 'Location', 'northeast');
set(gca, 'FontSize', fS, 'FontName', 'Arial', 'XTick', (0 : 10 : 40), ...
                                                    'YTick', (0 : 5 : 20)); 
n10 = subplot(5, 2, 10);
hold on 
for i = 1 : numberOfPhases 
    plot(G.cells.centroids(:, 1), 100 * (1 - fluid.K(porosity - ...
     states{phaseEndSteps(i)}.b - states{phaseEndSteps(i)}.c) ...
        ./ K0), 'color', pltCls{i}, 'LineWidth', lW, ...
                                                       'LineStyle', lS{i});                   
end
for i = 1 : numberOfPhases 
    k = 100 * (1 - fluid.K(porosity - states_OPM{i + 1}.b - ...
                                               states_OPM{i + 1}.c) ./ K0);  
    plot(G.cells.centroids(1 : 2 : end, 1), k(1 : 2 : end), 'o', ...
        'MarkerSize', mS, 'MarkerEdgeColor', [0 0 1], 'LineStyle', 'none');
end
line([pL(1) pL(1)], [0 100], 'color', 'red', 'LineStyle', ':', ...
                                                          'LineWidth', lW);                                                  
line([pL(2) pL(2)], [0 100], 'color', 'red', 'LineStyle', ':', ...
                                                          'LineWidth', lW);
line([pL(1) pL(2)], [100 100], 'color', 'red', 'LineStyle', ':', ...
                                                          'LineWidth', lW);
line([pL(1) pL(2)], [0 0], 'color', 'red', 'LineStyle', ':', ...
                                                          'LineWidth', lW);
hold off
xlim([0 40]);
ylim([0 100]);
xlabel('x [m]', 'FontSize', fS, 'Interpreter', 'latex');        
ylabel('$|\Delta K/K_0|$ [$\%$]', 'FontSize', fS, 'Interpreter', 'latex');
ax = gca;
ax.YAxis.Exponent = 0;
grid on
title('(g) Permeability', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
cb = legend('$20\;$h', '$40\;$h', ...
                           '$140\;$h$\qquad \qquad \qquad$','$160\;$h', ...
                           '$180\;$h','$230\;$h$\qquad \qquad \qquad$', ...
                                    '$250\;$h', '$270\;$h', '$500\;$h', ...
                                     'Location', 'best', 'Interpreter', ...
                                                  'latex', 'FontSize', fS);
set(cb, 'position', [0.1 0.78 0.8 0.15]);
cb.NumColumns = 3;
set(gca,'FontSize', fS, 'FontName', 'Arial', 'XTick', (0 : 10 : 40), ...
                                                  'YTick', (0 : 20 : 100));
ax8 = axes('position', [0.89 0.225 0.01 0.01], 'XColor', 'none', ...
                                                          'YColor','none');
hold on
plot(0, 0, 'color', [0 0 0], 'LineStyle', '-'); 
plot(-10, -200, 'o', 'MarkerSize', mS, 'MarkerEdgeColor', [0 0 1], ...
                                                      'LineStyle', 'none');
plot(-10, -200, 'o', 'MarkerSize', mS, 'MarkerEdgeColor', [1 1 1], ...
                                                      'LineStyle', 'none');                                                  
hold off
legend('MRST', 'OPM', 'Location', 'northeast');
set(gca, 'FontSize', fS, 'FontName', 'Arial', 'XTick', (0 : 10 : 40), ...
                                                  'YTick', (0 : 20 : 100));
%print -depsc2 Fig3.eps    