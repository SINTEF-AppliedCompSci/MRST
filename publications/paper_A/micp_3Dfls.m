% Setting up and solving the 3D flow leaky system (3Dfls).
% In MATLAB, this file produces Figures 12 and 13 in:
%
% Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E., 2021. Practical 
% approaches to study microbially induced calcite precipitation at the 
% field scale. Int. J. Greenh. Gas Control 106, 103256.
% https://doi.org/10.1016/j.ijggc.2021.103256. 
%
% In GNU Octave, this file creates and prints the results in the folder 
% vtk_micp_3Dfls which can be visualized using ParaView. The example 
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

% To get distmesh for first time run the following lines
% pth = fullfile(ROOTDIR,'utils','3rdparty','distmesh');
% mkdir(pth)
% unzip('http://persson.berkeley.edu/distmesh/distmesh.zip', pth);
% mrstPath('reregister','distmesh', pth);

% Required modules
pth = fullfile(ROOTDIR, 'utils', '3rdparty', 'distmesh');
mrstPath('reregister', 'distmesh', pth);
mrstModule add ad-blackoil ad-core ad-micp distmesh

% Grid
L = 500;  % Reservoir half length/width, m
H = 160;  % Reservoir heigth, m
nz = 0.5;
if exist('OCTAVE_VERSION', 'builtin') ~= 0 % GNU Octave
    hmin = 0.3;
    hmid = 15;
    hmax = 500;
    fd = @(p) drectangle(p, -L, L, -L, L);
    fh = @(p) min(min(1 + 0.58 * abs(dcircle(p, -100, 0, 0)), hmid) .* ...
     (abs(dcircle(p, -100, 0, 0)) < 50) + min(hmid + 0.6 * abs(dcircle( ...
        p, -100, 0, 50)), hmax) .* (abs(dcircle(p, -100, 0, 0)) >= 50), ...
              min(hmin + 0.6 * abs(dcircle(p, 0, 0, 0)), hmid) .* (abs( ...
                     dcircle(p, 0, 0, 0)) < 50) + min(hmid + 0.6 * abs( ...
                          dcircle(p, 0, 0, 50)), hmax) .* (abs(dcircle( ...
                                                      p, 0, 0, 0)) >= 50));
    [p, t] = distmesh2d(fd, fh, hmin, ...
                 [-L, -L ; L , L], [-L, -L ; L, -L ; -L, L ; L, L ; 0, 0]);
    close
    G = makeLayeredGrid(pebi(triangleGrid(p, t)), nz * H);
else % MATLAB
    hmax = 30;
    fd = @(p) drectangle(p, -L, L, -L, L);
    [p, t] = distmesh2d(fd, @huniform, hmax, [-L, -L ; L, L], ...
                                          [-L, -L ; L, -L ; -L, L ; L, L]);
    close
    Pw = [];
    for l = 50 * exp(-3 : 0.125 : 0)
        [x, y, z] = cylinder(l, 28); 
        Pw = [Pw [x(1,:); y(1,:)]];
    end
    Pw = [Pw [0 ; 0]];
    Pw1 = bsxfun(@plus, Pw, [-100 ; 0]);
    Pw = [];
    for l = 50 * exp(-5.1 : 0.125 : 0)
        [x, y, z] = cylinder(l, 28); 
        Pw = [Pw [x(1, :); y(1, :)]];
    end
    Pw = [Pw [0 ; 0]];
    Pw2 = bsxfun(@plus, Pw, [0 ; 0]);
    P = unique([Pw1' ; Pw2' ; p(:,1) p(:,2) ; 0 0], 'rows');
    G = makeLayeredGrid(pebi(triangleGrid(P)),nz * H);
end
rf= -4;
rs = 4 / (nz * 30);
rr = rf : rs : 0;
rfu= -4;
rsu = 4 / (nz * 30);
rru = rfu : rsu : 0;
mm = G.nodes.num / (nz * H + 1);
h1 = 30 * exp(rru);
for i = 0 : nz * 30
    G.nodes.coords(1 + mm * i : 1 : mm * (1 + i), 3) = ...
                                                   ones(mm, 1) * h1(i + 1);
end
for i = nz * 30 + 1 : nz * 130 - 1
    G.nodes.coords(1 + mm * i : mm * (1 + i), 3) = ...
                         G.nodes.coords(1 + mm * i : mm * (1 + i), 3) / nz;
end
for i = nz * 130 : nz * 160
    G.nodes.coords(1 + mm * i : mm * (1 + i), 3) = ...
            (G.nodes.coords(1 + mm * i : mm * (1 + i), 3) - nz * 130) * ...
                                       exp(rr(1 + i - nz * 130)) /nz + 130;
end
G = computeGeometry(G);
c = G.cells.centroids;
G = removeCells(G, (c(:, 1) < -50 * exp(-5.1) / 2 | c(:, 1) > ...
                     50 * exp(-5.1) / 2) & (c(:, 3) < 130 & c(:, 3) > 30));
c = G.cells.centroids;
G = removeCells(G, (c(:, 3) < 130 & c(:, 3) > 30) & (c(:, 2) <- 0.1 | ...
                                                           c(:, 2) > 0.1));
G = computeGeometry(G);
c = G.cells.centroids;
C = ones(G.cells.num, 1);

% Rock
K0 = 2e-14 * C;              % Aquifer permeability, m^2
cellsfrac =  G.cells.indexMap;
cellsfrac1 = cellsfrac((c(:, 1) > -50 * exp(-5.1) / 2 & c(:, 1) < 50 * ...
             exp(-5.1) / 2) & (c(:, 2) < 50 * exp(-5.1) / 2 & c(:, 2) > ...
                                                     -50 * exp(-5.1) / 2));
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

M(1, :)  = [15 * hour,   dt,    3   , 0.01,   0,      0];
M(2, :)  = [11 * hour,   dt,    3   , 0,      0,      0];
M(3, :)  = [74 * hour,   dt,    0   , 0,      0,      0];
M(4, :)  = [30 * hour,   dt,    3   , 0,      0.04,   0];
M(5, :)  = [5 * hour,    dt,    3   , 0,      0,      0];
M(6, :)  = [25 * hour,   dt,    0   , 0,      0,      0];
M(7, :)  = [40 * hour,   dt,    3   , 0,      0,    300];
M(8, :)  = [10 * hour,   dt,    3   , 0,      0,      0];
M(9, :)  = [390 * hour,  dt,    0   , 0,      0,      0];
M(10, :) = [30 * hour,   dt,    3   , 0,      0.04,   0];
M(11, :) = [20 * hour,   dt,    3   , 0,      0,      0];
M(12, :) = [20 * hour,   dt,    0   , 0,      0,      0];
M(13, :) = [20 * hour,   dt,    3   , 0,      0,    300];
M(14, :) = [20 * hour,   dt,    3   , 0,      0,      0];
M(15, :) = [90 * hour,   dt,    0   , 0,      0,      0];
M(16, :) = [20 * hour,   dt,    3   , 0,      0,    300];
M(17, :) = [20 * hour,   dt,    3   , 0,      0,      0];
M(18, :) = [110 * hour,  dt,    0   , 0,      0,      0];                  

% Create Well
r = 0.15;
Whu = 1 / 10;    
Whb = 1 - Whu;
[~, iw] = min(abs((c(:, 1) + 100) .^ 2 + c(:, 2) .^ 2));
cellsW = 1 : G.cells.num;
cellsWu = cellsW(abs(c(:, 1) - c(iw, 1)) < 0.01 & abs(c(:, 2) - ...
                          c(iw, 2)) < 0.01 & c(:,3) > 130 & c(:, 3) < 133);
W = addWell([],G, rock, cellsWu, 'Type', 'rate', 'Comp_i', [1, 0], ...
                                        'Val', Whu * M(1, 3), 'Radius', r);
cellsWb = cellsW(abs(c(:, 1) - c(iw, 1)) < 0.01 & abs(c(:, 2) - ...
                                         c(iw, 2)) < 0.01 & c(:, 3) > 133);
W = addWell(W, G, rock, cellsWb, 'Type', 'rate', 'Comp_i', [1, 0], ...
                                        'Val', Whb * M(1, 3), 'Radius', r);
for i = 1 : 2
    W(i).o = 0;
    W(i).u = 0;
    W(i).m = 0;
end
W(1).m = M(1, 4); 
W(1).o = M(1, 5);  
W(1).u = M(1, 6); 
G.injectionwellonboundary = 0;
             
% Gravity
gravity on

% Boundary Condition
f = boundaryFaces(G);
f = f(abs(G.faces.normals(f, 1)) > eps & (G.faces.centroids(f, 1) < ...
                                -L + 2 | G.faces.centroids(f, 1) > L - 2));
fp = G.faces.centroids(f, 3) * fluid.rhoWS * norm(gravity);
bc = addBC([], f, 'pressure', fp, 'sat', [0 0]);
bc.m = zeros(size(bc.sat, 1), 1);
bc.o = zeros(size(bc.sat, 1), 1);
bc.u = zeros(size(bc.sat, 1), 1);
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

% Initial Condition
state0 = initState(G, W, c(:, 3) * fluid.rhoWS * norm(gravity), [1, 0]);
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
    fn = getPlotAfterStepMICP(state0, model, 340, 20);
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
QCO2 = 1600 / day;
cellsW = 1 : G.cells.num;
cellsW = cellsW(abs(c(:, 1) - c(iw, 1)) < 0.01 & abs(c(:, 2) - ...
                                         c(iw, 2)) < 0.01 & c(:, 3) > 130);
W = addWell([], G, rock, cellsW, 'Type', 'rate', 'Comp_i', [0, 1], ... 
                                                 'Val', QCO2, 'Radius', r);
   
% Make schedule
schedule_co2 = simpleSchedule(timesteps, 'W', W, 'bc', bc); 

% Initial state
state0 = initState(G, W, c(:, 3) * fluid.rhoWS * norm(gravity), [1, 0]);

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
cellsfa =  1 : G.faces.num;
cellsfac = cellsfa(G.faces.centroids(:, 3) < 80 + 1000 * eps & ...
                                G.faces.centroids(:, 3) > 80 - 1000 * eps);
if exist('OCTAVE_VERSION', 'builtin') == 0
    fn = getPlotAfterStepCO2(state0, model, 340, 20);
    [~, statesco2] = simulateScheduleAD(state0, modela, schedule_co2, ...
                                                        'afterStepFn', fn);
    for i = 1 : ntco2
        lr0(i) = abs(statesco2{i}.flux(cellsfac(1), 2));
    end 
    statese = statesco2{end};
    clear statesco2
    [~, statesco2] = simulateScheduleAD(state0, modelb, schedule_co2, ...
                                                        'afterStepFn', fn);
    for i = 1 : ntco2
        lr1(i) = abs(statesco2{i}.flux(cellsfac(1), 2));
    end 
    statesf = statesco2{end};
    clear statesco2                                                 
    [~, statesco2] = simulateScheduleAD(state0, modelc, schedule_co2, ...
                                                        'afterStepFn', fn);
    for i = 1 : ntco2
        lr2(i) = abs(statesco2{i}.flux(cellsfac(1), 2));
    end 
    statesg = statesco2{end};
    clear statesco2                                                
    [~, statesco2] = simulateScheduleAD(state0, modeld, schedule_co2, ...
                                                        'afterStepFn', fn); 
    for i = 1:ntco2
        lr3(i) = abs(statesco2{i}.flux(cellsfac(1), 2));
    end 
    statesh = statesco2{end};
    clear statesco2                                                 
else
    [~, statesco2] = simulateScheduleAD(state0, modela, schedule_co2);
    for i = 1 : ntco2
        lr0(i) = abs(statesco2{i}.flux(cellsfac(1), 2));
    end 
    statese=statesco2{end};
    clear statesco2
    [~, statesco2] = simulateScheduleAD(state0, modelb, schedule_co2);
    for i = 1 : ntco2
        lr1(i) = abs(statesco2{i}.flux(cellsfac(1), 2));
    end 
    statesf = statesco2{end};
    clear statesco2 
    [~, statesco2] = simulateScheduleAD(state0, modelc, schedule_co2);
    for i = 1 : ntco2
        lr2(i) = abs(statesco2{i}.flux(cellsfac(1), 2));
    end 
    statesg = statesco2{end};
    clear statesco2
    [~, statesco2] = simulateScheduleAD(state0, modeld, schedule_co2);
    for i = 1 : ntco2
        lr3(i) = abs(statesco2{i}.flux(cellsfac(1), 2));
    end 
    statesh = statesco2{end};
    clear statesco2 
end

% Write the results to be read in ParaView (GNU Octave)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mkdir vtk_micp_3Dfls;
    cd vtk_micp_3Dfls;
    mrsttovtk(G, states, 'states_3Dfls', '%f');
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
    Ki = fluid.K(porosity - states{i}.c - states{i}.b);
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