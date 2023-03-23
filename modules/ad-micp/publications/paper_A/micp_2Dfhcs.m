% Setting up and solving the 2D flow horizontal circular system (2Dfhcs).
% In MATLAB, this file produces Figure 6 in:
%
% Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E., 2021. Practical 
% approaches to study microbially induced calcite precipitation at the 
% field scale. Int. J. Greenh. Gas Control 106, 103256.
% https://doi.org/10.1016/j.ijggc.2021.103256. 
%
% In GNU Octave, this file creates and prints the results in the folder 
% vtk_micp_2Dfhcs which can be visualized using ParaView. The example 
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

% To get distmesh for first time, uncomment and run the following lines
% pth = fullfile(ROOTDIR,'utils','3rdparty','distmesh');
% mkdir(pth)
% unzip('http://persson.berkeley.edu/distmesh/distmesh.zip', pth);
% mrstPath('reregister','distmesh', pth);

% Required modules
pth = fullfile(ROOTDIR, 'utils', '3rdparty', 'distmesh');
mrstPath('reregister', 'distmesh', pth);
mrstModule add ad-blackoil ad-core ad-micp distmesh
                                                               
% Grid 
R = 75;        % Reservoir radius, m  
B = 25;        % hmin to hmax transition radius, m
hmin = 0.75;   % Minimum grid size, m
hmax = 10;     % Maximum grid size, m
fd = @(p) dcircle(p, 0, 0, R);
fh = @(p) min(hmin + 0.3 * abs(dcircle(p, 0, 0, 0)), hmin) .* ...
 (abs(dcircle(p, 0, 0, 0)) < B) + min(hmin + 0.3 * abs(dcircle(p, 0, 0, ... 
                             B)), hmax) .* (abs(dcircle(p, 0, 0, 0)) >= B);
[p, t] = distmesh2d(fd, fh, hmin, [-R, -R ; R , R], ...
                                [-R , -R ; R , -R ; -R , R ; R, R ; 0, 0]);
G = makeLayeredGrid(pebi(triangleGrid(p, t)), 1);
close
G = computeGeometry(G);
c = G.cells.centroids;
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
dt_on = hour; % Time step when the well is on
dt_off = 10 * hour;  % Time step when the well is off

M(1, :) = [20 * hour,   dt_on,  1.2e-3  , 0.01, 0,      0];
M(2, :) = [20 * hour,   dt_on,  1.2e-3  , 0,    0,      0];
M(3, :) = [100 * hour,  dt_off, 0       , 0,    0,      0];
M(4, :) = [20 * hour,   dt_on,  1.2e-3  , 0,    0.04,   0];
M(5, :) = [20 * hour,   dt_on,  1.2e-3  , 0,    0,      0];
M(6, :) = [50 * hour,   dt_off, 0       , 0,    0,      0];
M(7, :) = [20 * hour,   dt_on,  1.2e-3  , 0,    0,    300];
M(8, :) = [20 * hour,   dt_on,  1.2e-3  , 0,    0,      0];
M(9, :) = [230 * hour,  dt_off, 0       , 0,    0,      0];              

% Create Well
r = 0.15;
cellsW =  G.cells.indexMap;
cellsW = cellsW(c(:, 1) .^ 2 + c(:, 2) .^ 2 < (hmin / 2) ^ 2);
W = addWell([], G, rock, cellsW, 'Type', 'rate', 'Comp_i', [1, 0], ...
                                              'Val', M(1, 3), 'Radius', r);
G.injectionwellonboundary = 0; 
W.m = M(1, 4);      
W.o = M(1, 5);
W.u = M(1, 6);                                               

% Boundary condition
f = boundaryFaces(G);
f = f((abs(G.faces.normals(f, 2)) > eps | abs(G.faces.normals(f, 1)) > ...
   eps) & G.faces.centroids(f, 1) .^ 2 + G.faces.centroids(f, 2) .^ 2 > ...
                                                             R - hmax / 4);
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
    mkdir vtk_micp_2Dfhcs;
    cd vtk_micp_2Dfhcs;
    mrsttovtk(G, states, 'states_2Dfhcs', '%f');
    return
end

% Figure 6 paper (MATLAB)
fS = 8;
figure;
c = flipud(jet);
sz = size(c, 1);
c = c((round(70 * sz / 256)) : end, :);
cc(:, 1) = (0.75 : 0.01 : 1)';
cc(:, 2) = (0.75 : 0.01 : 1)';
cc(:, 3) = (0.75 : -0.03 : 0)';
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6.83 4]);
n1 = subplot(1, 2, 1);
colormap (n1, cc);
caxis([0 1e-12]);
cb = colorbar; 
title(cb, '$m^2$', 'FontSize', fS, 'Interpreter', 'latex', ...
                                                      'FontName', 'Arial');
set(cb, 'position', [0.475 0.25 0.01 0.55], 'Ticks', [0 1e-12], ...
                     'FontSize', fS, 'TickLabels', {'0', '$10^{-12}$'}, ...
                                          'TickLabelInterpreter', 'latex');
plotCellData(G, K0);
title('Initial permeability', 'Interpreter', 'latex', 'FontSize', fS);
axis equal tight;
xlabel({'x [m]', '{(a)}'});
ylabel('y [m]');
xlim([-R, R]);
ylim([-R, R]);
view(0, 90)
set(gca, 'FontSize', fS, 'XTick', -R : 25 : R, 'YTick', ...
                        -R : 25 : R, 'color', 'none', 'FontName', 'Arial');
n2 = subplot(1, 2, 2);
s = plotCellData(G, 100 * (K0 - fluid.K(porosity - states{end}.c - ...
                                                    states{end}.b)) ./ K0);
title('Permeability (phase I MICP)', 'Interpreter', 'latex', ...
                                                           'FontSize', fS);
axis equal tight;
colormap(n2, c)
caxis([0 100]);
cb = colorbar; 
title(cb, '$\%$', 'Interpreter', 'latex');
set(cb,'position', [0.925 0.25 0.01 0.55], 'YTick', 0 : 20 : 100, ...
                                                           'FontSize', fS);
xlabel({'x [m]', '(b)'});
ylabel('y [m]');
s.EdgeColor = 'none';
xlim([-R, R]);
ylim([-R, R]);
view(0, 90)
grid off;
set(gca, 'FontSize', fS, 'XTick', -R : 25 : R, 'YTick', ...
                        -R : 25 : R, 'color', 'none', 'FontName', 'Arial');
rectangle('Position', [10, -2.5 , 5, 5], 'LineWidth', 2 , ...
                                 'LineStyle', '-', 'edgecolor', '[0 0 0]');
ax3 = axes('position', [0.15 0.78 0.05 0.06]);
box on 
axis([-1 1 -1 1]);
xlim([-1 1])
ylim([-1 1])
colormap(ax3, 'parula')
s = plotGrid(G);
view(0, 90)
set(gca, 'FontSize', fS, 'XTick', -1 : 1, 'YTick', -1 : 1, ...
                                     'color', 'none', 'FontName', 'Arial');
%print -depsc2 Fig6.eps