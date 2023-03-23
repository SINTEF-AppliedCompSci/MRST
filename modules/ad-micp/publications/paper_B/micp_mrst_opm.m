% Setting up and solving the system to compare to OPM simulation results.
%
% In MATLAB, this file produces Figure 3 in:
%
% Landa-Marbán, D., Kumar, K., Tveit, S., Gasda, S.E. Numerical studies of 
% CO2 leakage remediation by micp-based plugging technology. In: Røkke, 
% N.A. and Knuutila, H.K. (Eds) Short Papers from the 11th International 
% Trondheim CCS conference, ISBN: 978-82-536-1714-5, 284-290. 
%
% In GNU Octave, this file creates and prints the results in the folder 
% vtk_micp_mrst_opm which can be visualized using ParaView. The example 
% assumes MRST is the Matlab/Octave path. For information on 
% MRST-functions, confer the MRST documentation at:
%
% http://www.sintef.no/projectweb/mrst/
%
%{
Copyright 2021, NORCE Norwegian Research Centre AS, Computational 
Geosciences and Modeling. 

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}

% Required modules
mrstModule add ad-blackoil ad-core ad-micp

% Grid 
L = 100;                     % Aquifer length, m
G = tensorGrid(0 : L, [0 1], [0 1]);
G = computeGeometry(G);
C = ones(G.cells.num, 1); 

% Rock
K0 = 1e-14 * C;              % Aquifer permeability, m^2
porosity = 0.15;             % Aquifer porosity, [-]
rock = makeRock(G, K0, porosity);

% Fluid properties
fluid.muw = 2.535e-4;        % Water viscocity, Pa s                            
fluid.bW = @(p) 0 * p + 1;   % Water formation volume factor, [-]
fluid.bO = @(p) 0 * p + 1;   % CO2 formation volume factor, [-]
fluid.rhoWS = 1045;          % Water density, kg/m^3
fluid.rhoOS = 479;           % CO2 density, kg/m^3

% Remaining model parameters (we put them on the fluid structure)
fluid.rho_b = 35;            % Density (biofilm), kg/m^3
fluid.rho_c = 2710;          % Density (calcite), kg/m^3
fluid.k_str = 2.6e-10;       % Detachment rate, m/(Pa s)
fluid.diffm = 0;             % Diffusion coefficient (microbes), m^2/s
fluid.diffo = 0;             % Diffusion coefficient (oxygen), m^2/s
fluid.diffu = 0;             % Diffusion coefficient (urea), m^2/s
fluid.alphaL = 0;            % Disperison coefficient (longitudinal), m
fluid.alphaT = 0;            % Disperison coefficient (transverse), m
fluid.eta = 3;               % Fitting factor, [-]
fluid.k_o = 2e-5;            % Half-velocity constant (oxygen), kg/m^3
fluid.k_u = 21.3;            % Half-velocity constant (urea), kg/m^3
fluid.mu = 4.17e-5;          % Maximum specific growth rate, 1/s
fluid.mu_u = 0.0161;         % Maximum rate of urease utilization, 1/s
fluid.k_a = 8.51e-7;         % Microbial attachment rate, 1/s                                         
fluid.k_d = 3.18e-7;         % Microbial death rate, 1/s
fluid.Y = 0.5;               % Yield growth coefficient, [-]
fluid.Yuc = 1.67;            % Yield coeccifient (calcite/urea), [-]
fluid.F = 0.5;               % Oxygen consumption factor, [-]
fluid.crit = 0.1;            % Critical porosity, [-]
fluid.kmin = 1e-20;          % Minimum permeability, m^2
fluid.cells = C;             % Array with all cells, [-]
fluid.ptol = 1e-4;           % Porosity tolerance to stop the simulation  

% Porosity-permeability relationship
fluid.K = @(poro) (K0 .* ((poro - fluid.crit) / (porosity - fluid.crit))...
               .^ fluid.eta + fluid.kmin) .* K0 ./ (K0 + fluid.kmin) .* ...
                  (poro > fluid.crit) + fluid.kmin .* (poro <= fluid.crit);      
              
% Injection strategy
N = 9;              % Number of injection phases in the injection strategy.
M = zeros(N, 6);    % The entries per row are: time, dt, rate, m, o, and u.
dt_on = hour;       % Time step when the well is on.
dt_off = 5 * hour;  % Time step when the well is off.

M(1, :) = [20 * hour,   dt_on,  2 / day , 0.01, 0,      0];
M(2, :) = [20 * hour,   dt_on,  2 / day , 0,    0,      0];
M(3, :) = [100 * hour,  dt_off, 0       , 0,    0,      0];
M(4, :) = [20 * hour,   dt_on,  2 / day , 0,    0.04,   0];
M(5, :) = [20 * hour,   dt_on,  2 / day , 0,    0,      0];
M(6, :) = [50 * hour,   dt_off, 0       , 0,    0,      0];
M(7, :) = [20 * hour,   dt_on,  2 / day , 0,    0,    300];
M(8, :) = [20 * hour,   dt_on,  2 / day , 0,    0,      0];
M(9, :) = [230 * hour,  dt_off, 0       , 0,    0,      0];

% Create well
r = 0.15;                       % Well radius, m
W = addWell([], G, rock, 1, 'Type', 'rate', 'Comp_i', [1, 0], ...
                                              'Val', M(1, 3), 'Radius', r);
W = addWell(W, G, rock, G.cells.num, 'Type', 'bhp', 'Comp_i', [1, 0], ...
                                                  'Val', atm, 'Radius', r);

G.injectionwellonboundary = 1;  % The injection well is on the boundary
G.cellsinjectionwell = 1; 
                                              
for i = 1 : 2
    W(i).m = 0;
    W(i).o = 0;
    W(i).u = 0;
end
W(1).m = M(1, 4);      
W(1).o = M(1, 5);
W(1).u = M(1, 6);

% Setup some schedule
nt = sum(M(:, 1) ./ M(:, 2)); 
timesteps = repmat(dt_on, nt, 1);
schedule = simpleSchedule(timesteps, 'W', W);
for i = 2 : N 
    schedule.control(i) = schedule.control(i - 1);
    schedule.step.control(sum(M(1 : i - 1, 1) ./ M(1 : i - 1, 2)) + 1 : ...
                                                                  end) = i;
    schedule.step.val(sum(M(1 : i - 1, 1) ./ M(1 : i - 1, 2)) + 1 : ...
                                                            end) = M(i, 2);
    schedule.control(i).W(1).val = M(i, 3);
    schedule.control(i).W(1).m = M(i, 4);
    schedule.control(i).W(1).o = M(i, 5);
    schedule.control(i).W(1).u = M(i, 6);
end

% Maximum injected oxygen and urea concentrations.
fluid.Comax = max(M(:, 5));             
fluid.Cumax = max(M(:, 6));

% Create model
model = MICPModel(G, rock, fluid);
model.toleranceMB = 1e-14;
model.nonlinearTolerance = 1e-14;

% Set up solver
solver = getNonLinearSolver(model);
solver.LinearSolver.tolerance = 1e-14;

% Initial condition
state0 = initState(G, W, atm, [1, 0]);
state0.m = zeros(G.cells.num, 1);
state0.o = zeros(G.cells.num, 1);
state0.u = zeros(G.cells.num, 1);
state0.b = zeros(G.cells.num, 1);
state0.c = zeros(G.cells.num, 1);

% Simulate case (GNU Octave/MATLAB)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    ok = 'true';
    fn = checkCloggingMICP(ok);
else
    fn = getPlotAfterStepMICP(state0, model, 0, 270);
end
[~, states] = simulateScheduleAD(state0, model, schedule, ...
                             'NonLinearSolver', solver, 'afterStepFn', fn);

% Write the results to be read in ParaView (GNU Octave)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mkdir vtk_micp_mrst_opm;
    cd vtk_micp_mrst_opm;
    mrsttovtk(G, states, 'states', '%f');
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
for i = 1 : N 
    plot(G.cells.centroids(1 : 2 : end, 1), ...
            states_OPM{i + 1}.m(1 : 2 : end, 1), 'o', 'MarkerSize', mS, ...
                          'MarkerEdgeColor', [0 0 1], 'LineStyle', 'none');
    plot(G.cells.centroids(:, 1), ...
         states{sum(M(1 : i, 1) ./ M(1 : i, 2))}.m, 'color', pltCls{i}, ...
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
for i = 1 : N 
    plot(G.cells.centroids(1 : 2 : end, 1), ...
            states_OPM{i + 1}.o(1 : 2 : end, 1), 'o', 'MarkerSize', mS, ...
                          'MarkerEdgeColor', [0 0 1], 'LineStyle', 'none');
    plot(G.cells.centroids(:, 1), ...
         states{sum(M(1 : i, 1) ./ M(1 : i, 2))}.o, 'color', pltCls{i}, ...
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
for i = 1 : N 
    plot(G.cells.centroids(1 : 2 : end, 1), ...
            states_OPM{i + 1}.b(1 : 2 : end, 1), 'o', 'MarkerSize', mS, ...
                          'MarkerEdgeColor', [0 0 1], 'LineStyle', 'none');
    plot(G.cells.centroids(:, 1), ...
         states{sum(M(1 : i, 1) ./ M(1 : i, 2))}.b, 'color', pltCls{i}, ...
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
for i = 1 : N 
    plot(G.cells.centroids(1 : 2 : end, 1), ...
            states_OPM{i + 1}.u(1 : 2 : end, 1), 'o', 'MarkerSize', mS, ...
                          'MarkerEdgeColor', [0 0 1], 'LineStyle', 'none');
    plot(G.cells.centroids(:, 1), ...
         states{sum(M(1 : i, 1) ./ M(1 : i, 2))}.u, 'color', pltCls{i}, ...
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
for i = 1 : N 
    plot(G.cells.centroids(1 : 2 : end, 1), ...
            states_OPM{i + 1}.c(1 : 2 : end, 1), 'o', 'MarkerSize', mS, ...
                          'MarkerEdgeColor', [0 0 1], 'LineStyle', 'none');
    plot(G.cells.centroids(:, 1), ...
         states{sum(M(1 : i, 1) ./ M(1 : i, 2))}.c, 'color', pltCls{i}, ...
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
for i = 1 : N 
    plot(G.cells.centroids(1 : 2 : end, 1), 100 / porosity * ( ...
              states_OPM{i + 1}.b(1 : 2 : end, 1) + states_OPM{i + 1}.c ...
                              (1 : 2 : end, 1)), 'o', 'MarkerSize', mS, ...
                               'MarkerEdgeColor', [0 0 1], 'LineStyle', ...
                                                                   'none');
    plot(G.cells.centroids(:, 1), 100 / porosity * ( ...
                            states{sum(M(1 : i, 1) ./ M(1 : i, 2))}.b + ...
                            states{sum(M(1 : i, 1) ./ M(1 : i, 2))}.c), ... 
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
for i = 1 : N 
    plot(G.cells.centroids(:, 1), 100 * (1 - fluid.K(porosity - ...
     states{sum(M(1 : i, 1) ./ M(1 : i, 2))}.b - states{sum(M(1 : i, 1) ...
       ./ M(1 : i, 2))}.c) ./ K0), 'color', pltCls{i}, 'LineWidth', lW, ...
                                                       'LineStyle', lS{i});                   
end
for i = 1 : N 
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