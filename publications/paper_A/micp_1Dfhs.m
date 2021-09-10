% Setting up and solving the 1D flow horizontal system (1Dfhs).
% In MATLAB, this file produces Figure 5 in:
%
% Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E., 2021. Practical 
% approaches to study microbially induced calcite precipitation at the 
% field scale. Int. J. Greenh. Gas Control 106, 103256.
% https://doi.org/10.1016/j.ijggc.2021.103256. 
%
% In GNU Octave, this file creates and prints the results in the folder 
% vtk_micp_1Dfhs which can be visualized using ParaView. The example 
% assumes MRST is the Matlab/Octave path. For information on 
% MRST-functions, confer the MRST documentation at
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
L = 75;                  % Aquifer length, m
l = 25;                  % Region where MICP processes are more relevant, m
dw = 0.5;                % Size of the element for the well, m
dl = 0.05;               % Size of the elements inside l, m
X = [0 dw dw + dl : dl : l  L * exp(-1.075 : 0.025 : 0)];
G = tensorGrid(X, [0 1], [0 1]);
G = computeGeometry(G);
C = ones(G.cells.num, 1);

% Rock
K0 = 1e-12 * C;              % Aquifer permeability, m^2
porosity = 0.2;              % Aquifer porosity, [-]
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
fluid.diffm = 2.1e-9;        % Diffusion coefficient (microbes), m^2/s
fluid.diffo = 2.32e-9;       % Diffusion coefficient (oxygen), m^2/s
fluid.diffu = 1.38e-9;       % Diffusion coefficient (urea), m^2/s
fluid.alphaL = 1e-3;         % Disperison coefficient (longitudinal), m
fluid.alphaT = 4e-4;         % Disperison coefficient (transverse), m
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
N = 9; % Number of injection phases in the injection strategy
M = zeros(N, 6); % The entries per row are: time, dt, rate, m, o, and u.
dt_on = 20 * minute; % Time step when the well is on
dt_off = 5 * hour;   % Time step when the well is off

M(1, :) = [20 * hour,   dt_on,  2.4e-05 , 0.01, 0,      0];
M(2, :) = [20 * hour,   dt_on,  2.4e-05 , 0,    0,      0];
M(3, :) = [100 * hour,  dt_off, 0       , 0,    0,      0];
M(4, :) = [20 * hour,   dt_on,  2.4e-05 , 0,    0.04,   0];
M(5, :) = [20 * hour,   dt_on,  2.4e-05 , 0,    0,      0];
M(6, :) = [50 * hour,   dt_off, 0       , 0,    0,      0];
M(7, :) = [20 * hour,   dt_on,  2.4e-05 , 0,    0,    300];
M(8, :) = [20 * hour,   dt_on,  2.4e-05 , 0,    0,      0];
M(9, :) = [230 * hour,  dt_off, 0       , 0,    0,      0];
              
% Create well
r = 0.15;    
W = addWell([], G, rock, 1, 'Type', 'rate', 'Comp_i', [1, 0], ...
                                             'Val',  M(1, 3), 'Radius', r);
G.injectionwellonboundary = 1;
G.cellsinjectionwell = 1;
W.m = M(1, 4);      
W.o = M(1, 5);
W.u = M(1, 6); 

% Boundary condition
f = boundaryFaces(G);
f = f(abs(G.faces.normals(f, 1)) > 0 & G.faces.centroids(f, 1) > ...
                                                               X(end - 1));
bc = addBC([], f, 'pressure', atm, 'sat', [0 0]);
bc.m = zeros(size(bc.sat, 1), 1);
bc.o = zeros(size(bc.sat, 1), 1);
bc.u = zeros(size(bc.sat, 1), 1);
bc.b = zeros(size(bc.sat, 1), 1);
bc.c = zeros(size(bc.sat, 1), 1);

% Setup some schedule
nt = sum(M(:, 1) ./ M(:, 2));
timesteps = repmat(dt_on, nt, 1);
schedule = simpleSchedule(timesteps, 'W', W, 'bc', bc);
for i = 2 : N 
    schedule.control(i) = schedule.control(i - 1);
    schedule.step.control(sum(M(1 : i - 1, 1) ./ M(1 : i - 1, 2)) + 1 : ...
                                                                  end) = i;
    schedule.step.val(sum(M(1 : i - 1, 1) ./ M(1 : i - 1, 2)) + 1 : ...
                                                            end) = M(i, 2);
    schedule.control(i).W.val = M(i, 3);
    schedule.control(i).W.m = M(i, 4);
    schedule.control(i).W.o = M(i, 5);
    schedule.control(i).W.u = M(i, 6);
end

% Maximum injected oxygen and urea concentrations.
fluid.Comax = max(M(:, 5));             
fluid.Cumax = max(M(:, 6));  

% Create model
model = MICPModel(G, rock, fluid);
model.toleranceMB = 1e-14;
model.nonlinearTolerance = 1e-14;

% Solver
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
    fn = getPlotAfterStepMICP(state0, model, 0, 90);
end
[~, states] = simulateScheduleAD(state0, model, schedule, ...
                             'NonLinearSolver', solver, 'afterStepFn', fn);

% Write the results to be read in ParaView (GNU Octave)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mkdir vtk_micp_1Dfhs;
    cd vtk_micp_1Dfhs;
    mrsttovtk(G, states, 'states', '%f');
    return
end

% Figure 5 paper (MATLAB)
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
xlabel({'x [m]' ; '(a)'}, 'FontSize', fS, 'Interpreter', 'latex');        
ylabel('$c_m$ [kg/m$^3$]', 'FontSize', fS, 'Interpreter', 'latex');
grid on
title('Microbes', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca,'FontSize', fS, 'FontName', 'Arial', ...
                       'XTick', (0 : 10 : L), 'YTick', (0 : 0.002 : 0.01));
n2 = subplot(3, 3, 5);
hold on
for i = 1 : N 
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
xlim([0 L]);
ylim([0 0.04]);
xlabel({'x [m]' ; '(b)'}, 'FontSize', fS, 'Interpreter', 'latex');        
ylabel('$c_o$ [kg/m$^3$]', 'FontSize', fS, 'Interpreter', 'latex');
grid on
title('Oxygen', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca,'FontSize', fS, 'FontName', 'Arial', ...
                        'XTick', (0 : 10 : L), 'YTick', (0 : 0.01 : 0.04));
n3 = subplot(3, 3, 6);
hold on
for i = 1 : N 
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
    plot(G.cells.centroids(:, 1), ...
         states{sum(M(1 : i, 1) ./ M(1 : i, 2))}.b, 'color', pltCls{i}, ...
                                      'LineWidth', lW, 'LineStyle', lS{i});                   
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
xlabel({'x [m]' ; '(d)'}, 'FontSize', fS,'Interpreter', 'latex');        
ylabel('$\phi_b$ [$-$]', 'FontSize', fS, 'Interpreter', 'latex');
grid on
title('Biofilm', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', fS, 'FontName', 'Arial', ...
                    'XTick', (0 : 10 : L), 'YTick', (0 : 0.0001 : 0.0003));
n5 = subplot(3, 3, 8);
hold on
for i = 1 : N 
    plot(G.cells.centroids(:, 1), ...
         states{sum(M(1 : i, 1) ./ M(1 : i, 2))}.c, 'color', pltCls{i}, ...
                                      'LineWidth', lW, 'LineStyle', lS{i});                   
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
xlabel({'x [m]' ; '(e)'}', 'FontSize', fS, 'Interpreter', 'latex');        
ylabel('$\phi_c$ [$-$]', 'FontSize', fS, 'Interpreter', 'latex');
grid on
title('Calcite', 'FontSize', fS, 'FontName', 'Arial', ...
                                                   'Interpreter', 'latex');
set(gca, 'FontSize', fS,'FontName', 'Arial', ...
                         'XTick', (0 : 10 : L), 'YTick', (0 : 0.01 :0.05));                                                        
n6 = subplot(3, 3, 9);
hold on
for i = 1 : N 
    plot(G.cells.centroids(:, 1), 100 * (1 - fluid.K(porosity - ...
     states{sum(M(1 : i, 1) ./ M(1 : i, 2))}.b - states{sum(M(1 : i, 1) ...
       ./ M(1 : i, 2))}.c) ./ K0), 'color', pltCls{i}, 'LineWidth', lW, ...
                                                       'LineStyle', lS{i});                   
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
xlabel({'x [m]' ; '(f)'}', 'FontSize', fS, 'Interpreter', 'latex');        
ylabel('$|\Delta K/K_0|$ [$\%$]', 'FontSize', fS, 'Interpreter', 'latex');
grid on
title('Permeability', 'FontSize', fS, 'FontName', 'Arial', ...
                                                     'Interpreter','latex')
cb=legend('$t^I_1=\;\;20\;$h','$t^I_2=\;\;40\;$h', ...
               '$t^I_3=140\;$h$\qquad \qquad \qquad$','$t^I_4=160\;$h', ...
               '$t^I_5=180\;$h','$t^I_6=230\;$h$\qquad \qquad \qquad$', ...
                    '$t^I_7=250\;$h','$t^I_8=270\;$h','$t^I_9=500\;$h', ...
                   'Location','best','Interpreter','latex','FontSize', fS);
set(cb, 'position', [0.5 0.67 0.01 0.15]);
cb.NumColumns = 3;
set(gca, 'FontSize', fS, 'FontName', 'Arial', ...
                             'XTick',(0 : 10 : L), 'YTick',(0 : 20 : 100));
%print -depsc2 Fig5.eps