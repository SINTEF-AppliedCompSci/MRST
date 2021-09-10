% Setting up and solving the 2D flow leaky system (2Dfls).
% In MATLAB, this file produces Figures 10 and 11 in:
%
% Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E., 2021. Practical
% approaches to study microbially induced calcite precipitation at the
% field scale. Int. J. Greenh. Gas Control 106, 103256.
% https://doi.org/10.1016/j.ijggc.2021.103256.
%
% In GNU Octave, this file creates and prints the results in the folder
% vtk_micp_2Dfls which can be visualized using ParaView. The example
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
L = 500;        % Reservoir length, m
H = 160;        % Reservoir heigth, m
if exist('OCTAVE_VERSION', 'builtin') ~= 0 % GNU Octave
    Y = [0 : 0.25 : 135  H * exp(-0.16 : 0.02 : 0)];
    X = [0 1 50 * exp(-3.6 : 0.25 : 0) 100 - 50 * exp(0 : -0.25 : ...
       -4.75)  99.85 100.15 100.5 101 102 105 : 5 : 210 L * ...
                                                    exp(-0.85 : 0.05 : 0)];
    G = tensorGrid(X, Y, [0 1]);
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, (c(:, 1) < 100 - 0.3 | c(:, 1) > 100 + 0.3) & ...
                                           (c(:, 2) < 130 & c(:, 2) > 30));
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, (c(:, 1) < 99.9 | c(:, 1) > 100.1) & ...
                                           (c(:, 2) < 130 & c(:, 2) > 30));
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, c(:, 1) < -0.5 - eps & c(:, 2) > 130);
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, c(:, 1) < 99.8 & c(:, 2) < 30);
    G = computeGeometry(G);
else % MATLAB
    [X1, Y1] = meshgrid([-L : 10 : -50 180 :10 : L], [0 : 5 : 30 31 : ...
                                                         129 130 : 5 : H]);
    [X2, Y2] = meshgrid(-50 : 10 : 50, 0 : 5 : 30);
    [x, y] = meshgrid([100 - 0.3 100 100 + 0.3], 0 : 0.25 : H);
    [xc, yc] = meshgrid([-1 0 1], 130 : 0.25 : H);
    [xw, yw] = meshgrid([-50 * exp(0 : -0.25 : -3.6) 50 * ...
                                       exp(-3.6 : 0.25 : 0)], 130 : 1 : H);
    [xl, yl] = meshgrid([100 - 50 * exp(0 : -0.25 : -5) 100 + 80 * ...
                                              exp(-5 : 0.125 : 0)], 0 : H);
    [xwc1, ywc1] = meshgrid(50 * exp(-3.6 : 0.25 : 0), 130 : 0.25 : 137.5);
    [xwc2, ywc2] = meshgrid(50 * exp(-3.6 : 0.25 : 0), 137.5 : 0.5: 145);
    [xwc3, ywc3] = meshgrid(100 - 50 * exp(0 : -0.25 : -5), ...
                                                       130 : 0.25 : 137.5);
    [xwc4, ywc4] = meshgrid(100 - 50 * exp(0 : -0.25 : -5), ...
                                                        137.5 : 0.5 : 145);
    P = unique([X1(:) Y1(:); X2(:) Y2(:); x(:) y(:); xw(:) yw(:); xl(:) ...
          yl(:); xc(:) yc(:); xwc1(:) ywc1(:); xwc2(:) ywc2(:); xwc3(:) ...
                                        ywc3(:); xwc4(:) ywc4(:)], 'rows');
    G = triangleGrid(P);
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, (c(:, 1) < 100 - 0.3 | c(:, 1) > 100 + 0.3) & ...
                                           (c(:, 2) < 130 & c(:, 2) > 30));
    G = makeLayeredGrid(pebi(G), 1);
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, (c(:, 1) < 99.9 | c(:, 1) > 100.1) & ...
                                           (c(:, 2) < 130 & c(:, 2) > 30));
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, c(:, 1) < -0.5 - eps & c(:, 2) > 130);
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, c(:, 1) < 99.8 & c(:, 2) < 30);
    G = computeGeometry(G);
end
c = G.cells.centroids;
C = ones(G.cells.num,1);

% Rock
K0 = 2e-14 * C;              % Aquifer permeability, m^2
cellsfrac =  G.cells.indexMap;
cellsfrac1 = cellsfrac(c(:, 1) > 99.9 & c(:, 1) < 100.1 & ...
                                             c(:, 2) < 130 & c(:, 2) > 30);
cellsF =  G.cells.indexMap;
idx = ismember(cellsF, cellsfrac1);
K0(idx) = 1e-12;             % Leakage permeability, m^2
porosity = 0.15;             % Aquifer/leakage porosity, [-]
rock = makeRock(G, K0, porosity);

% Fluid properties
fluid.muw = 2.535e-4;        % Water viscocity, Pa s
fluid.muO = 3.95e-5;         % CO2 viscosity, Pa s
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
N = 18; % Number of injection phases in the injection strategy
M = zeros(N, 6); % The entries per row are: time, dt, rate, m, o, and u.
dt = hour;

M(1, :)  = [15 * hour,   dt,      6e-3  , 0.01,   0,      0];
M(2, :)  = [11 * hour,   dt,      6e-3  , 0,      0,      0];
M(3, :)  = [74 * hour,   dt,      0     , 0,      0,      0];
M(4, :)  = [30 * hour,   dt,      6e-3  , 0,      0.04,   0];
M(5, :)  = [5 * hour,    dt,      6e-3  , 0,      0,      0];
M(6, :)  = [25 * hour,   dt,      0     , 0,      0,      0];
M(7, :)  = [40 * hour,   dt,      6e-3  , 0,      0,    300];
M(8, :)  = [10 * hour,   dt,      6e-3  , 0,      0,      0];
M(9, :)  = [390 * hour,  dt,      0     , 0,      0,      0];
M(10, :) = [30 * hour,   dt,      6e-3  , 0,      0.04,   0];
M(11, :) = [20 * hour,   dt,      6e-3  , 0,      0,      0];
M(12, :) = [20 * hour,   dt,      0     , 0,      0,      0];
M(13, :) = [20 * hour,   dt,      6e-3  , 0,      0,    300];
M(14, :) = [20 * hour,   dt,      6e-3  , 0,      0,      0];
M(15, :) = [90 * hour,   dt,      0     , 0,      0,      0];
M(16, :) = [20 * hour,   dt,      6e-3  , 0,      0,    300];
M(17, :) = [20 * hour,   dt,      6e-3  , 0,      0,      0];
M(18, :) = [110 * hour,  dt,      0     , 0,      0,      0];

% Create well
r = 0.15;
Whu = 1 / 10;
Whb = 1 - Whu;
cellsW = 1 : G.cells.num;
cellsWu = cellsW(c(:, 1) < min(c(:, 1)) + 0.1 & c(:, 2) > 130 & ...
                                                            c(:, 2) < 133);
W = addWell([], G, rock, cellsWu, 'Type', 'rate', 'Comp_i', [1, 0], ...
                            'Val', Whu * M(1, 3), 'Radius', r, 'dir', 'y');
cellsWb = cellsW(c(:, 1) < min(c(:, 1)) + 0.1 & c(:, 2) > 133);
W = addWell(W, G, rock, cellsWb, 'Type', 'rate', 'Comp_i', [1, 0], ...
                            'Val', Whb * M(1, 3), 'Radius', r, 'dir', 'y');
for i = 1 : 2
    W(i).o = 0;
    W(i).u = 0;
    W(i).m = 0;
end
W(1).m = M(1, 4);
W(1).o = M(1, 5);
W(1).u = M(1, 6);
G.injectionwellonboundary = 1;
G.cellsinjectionwell = [cellsWu cellsWb];

% Gravity
gravity on
gravity y

% Boundary condition
f = boundaryFaces(G);
f = f(abs(G.faces.normals(f, 1)) > eps & (G.faces.centroids(f, 1) < -L ...
                                   + 2 | G.faces.centroids(f, 1) > L - 2));
fp = G.faces.centroids(f, 2) * fluid.rhoWS * norm(gravity);
bc = addBC([], f, 'pressure', fp, 'sat', [0 0]);
bc.o = zeros(size(bc.sat, 1), 1);
bc.u = zeros(size(bc.sat, 1), 1);
bc.m = zeros(size(bc.sat, 1), 1);
bc.b = zeros(size(bc.sat, 1), 1);
bc.c = zeros(size(bc.sat, 1), 1);

% Setup some schedule
nt = sum(M(:, 1) ./ M(:, 2));
timesteps = repmat(dt, nt, 1);
schedule = simpleSchedule(timesteps, 'W', W, 'bc', bc);
for i = 2 : N
    schedule.control(i) = schedule.control(i - 1);
    schedule.step.control(sum(M(1 : i - 1, 1) ./ M(1 : i - 1, 2)) + 1 : ...
                                                                  end) = i;
    schedule.step.val(sum(M(1 : i - 1, 1) ./ M(1 : i - 1, 2)) + 1 : ...
                                                            end) = M(i, 2);
    schedule.control(i).W(1).val = Whu * M(i, 3);
    schedule.control(i).W(2).val = Whb * M(i, 3);
    schedule.control(i).W(1).m = M(i, 4);
    schedule.control(i).W(1).o = M(i, 5);
    schedule.control(i).W(1).u = M(i, 6);
end

% Maximum injected oxygen and urea concentrations.
fluid.Comax = max(M(:, 5));
fluid.Cumax = max(M(:, 6));

% Create model
model = MICPModel(G, rock, fluid);

% Initial condition
state0 = initState(G, W, c(:, 2) * fluid.rhoWS * norm(gravity), [1, 0]);
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
                                                        'afterStepFn', fn);

% CO2 assesment
statesa = state0;
statesb = states{sum(M(1 : 9, 1) ./ M(1 : 9, 2))};
statesc = states{sum(M(1 : 15, 1) ./ M(1 : 15, 2))};
statesd = states{end};

% Setup some schedule
dt = day;
ntco2 = 100 * day / dt;
timesteps = repmat(dt, ntco2, 1);

% Create CO2 Well
QCO2 = (1600 / day) / L;
cellsW =  1 : G.cells.num;
cellsW = cellsW(c(:, 1) < min(c(:, 1)) + 0.1 & c(:,2) > 130);
W = addWell([], G, rock, cellsW, 'Type', 'rate', 'Comp_i', [0, 1], ...
                                     'Val', QCO2, 'Radius', r, 'dir', 'y');

% Make schedule
schedule_co2 = simpleSchedule(timesteps, 'W', W, 'bc', bc);

% Initial state
state0 = initState(G, W, G.cells.centroids(:, 2) * fluid.rhoWS * ...
                                                    norm(gravity), [1, 0]);

% Compute porosity and permeability
poro = porosity - statesa.c - statesa.b;
KK = fluid.K(poro);
rocka = makeRock(G, KK, poro);
poro = porosity - statesb.c - statesb.b;
KK = fluid.K(poro);
rockb = makeRock(G, KK, poro);
poro = porosity - statesc.c - statesc.b;
KK = fluid.K(poro);
rockc = makeRock(G, KK, poro);
poro = porosity - statesd.c - statesd.b;
KK = fluid.K(poro);
rockd = makeRock(G, KK, poro);

% Create model
modela = CO2Model(G, rocka, fluid);
modelb = CO2Model(G, rockb, fluid);
modelc = CO2Model(G, rockc, fluid);
modeld = CO2Model(G, rockd, fluid);

% Simulate
if exist('OCTAVE_VERSION', 'builtin') == 0
    fn = getPlotAfterStepCO2(state0, model, 0, 270);
    [~, statese] = simulateScheduleAD(state0, modela, schedule_co2, ...
                                                        'afterStepFn', fn);
    [~, statesf] = simulateScheduleAD(state0, modelb, schedule_co2, ...
                                                        'afterStepFn', fn);
    [~, statesg] = simulateScheduleAD(state0, modelc, schedule_co2, ...
                                                        'afterStepFn', fn);
    [~, statesh] = simulateScheduleAD(state0, modeld, schedule_co2, ...
                                                        'afterStepFn', fn);
else
    [~, statese] = simulateScheduleAD(state0, modela, schedule_co2);
    [~, statesf] = simulateScheduleAD(state0, modelb, schedule_co2);
    [~, statesg] = simulateScheduleAD(state0, modelc, schedule_co2);
    [~, statesh] = simulateScheduleAD(state0, modeld, schedule_co2);
end

% Compute leakage rate
cellsfa =  1 : G.faces.num;
cellsfac = cellsfa(G.faces.centroids(:, 2) < 80.6 & ...
         G.faces.centroids(:, 2) > 80.3 & abs(G.faces.normals(:, 2))> 0.1);
for i = 1 : ntco2
    lr0(i) = abs(statese{i}.flux(cellsfac(1), 2));
    lr1(i) = abs(statesf{i}.flux(cellsfac(1), 2));
    lr2(i) = abs(statesg{i}.flux(cellsfac(1), 2));
    lr3(i) = abs(statesh{i}.flux(cellsfac(1), 2));
end
statese = statese{end};
statesf = statesf{end};
statesg = statesg{end};
statesh = statesh{end};

% Write the results to be read in ParaView (GNU Octave)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mkdir vtk_micp_2Dfls;
    cd vtk_micp_2Dfls;
    mrsttovtk(G, states, 'states_2Dfls', '%f');
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
  Ki = fluid.K(porosity - states{i}.c - states{i}.b);
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
