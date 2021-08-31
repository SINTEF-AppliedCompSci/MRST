%% Workflow example for ad-micp
% This example aims to show complete workflow for creating, running, and
% analyzing a 2D flow system with a leakage path using the MICP
% mathematical model. To asses the co2 distribution before and after
% treatment, we use the TwoPhaseWaterGasModel in the MRST co2lab module (we
% use the system relationships/properties as in the basic_3D_example in
% the co2lab module).

% Required modules
mrstModule add deckformat ad-core ad-blackoil ad-micp ad-props mrst-gui ...
                                                                     co2lab

%% Reservoir geometry/properties and model parameters
% The domain has a length of 100 m and a height of 40 m. We remove the
% domain cells where there is caprock. The grid is coarser in order to run
% the example in couple of minutes. You should try with finer grids (see
% https://doi.org/10.1016/j.ijggc.2021.103256).

% Grid
L = 100;        % Reservoir length, m
ht = 4;         % Top aquifer heigth, m
hl = 32;        % Leakage heigth, m
hb = 4;         % Bottom aquifer heigth, m
H = ht+hl+hb;   % Reservoir heigth, m
D = 1500;       % Depth of aquifer top surface, m
a = 0.60;       % Leakage width, m
l = 15;         % Leakage distance from injection well, m
X = [0:l-1 l-a/2 l+a/2 l+1:21 L*exp(-1.5:.05:0)]; % X discretization
Z = [0:ht ht+1:1:hl+ht hl+ht+1:H];                % Z discretization
G = tensorGrid(X,[0 1],Z, 'depthz', repmat(D, 1, 2*size(X,2)));
G = computeGeometry(G);
c = G.cells.centroids;
G = removeCells(G,(c(:,1)<l-a/2 | c(:,1)>l+a/2) & c(:,3)<D+H-hb & ...
                                                              c(:,3)>D+ht);
G = computeGeometry(G);
c = G.cells.centroids;
C = ones(G.cells.num,1);

% Rock
K0 = 2e-14*C;                % Leakage permeability, m^2
cellsl = G.cells.indexMap;
cellsl = cellsl(c(:,3)<D+H-hb & c(:,3)>D+ht);
cellsF = G.cells.indexMap;
idx = ismember(cellsF,cellsl);
K0(idx) = 1e-12;             % Aquifer permeability, m^2
porosity = 0.15;             % Aquifer/leakage porosity, [-]
rock = makeRock(G, K0, porosity);

% Fluid properties
fluid.muw = 2.535e-4;        % Water viscocity, Pa s
fluid.bW = @(p) 0*p + 1;     % Water formation volume factor, [-]
fluid.rhoWS = 1045;          % Water density, kg/m^3

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
fluid.crit = .1;             % Critical porosity, [-]
fluid.kmin = 1e-20;          % Minimum permeability, m^2
fluid.cells = C;             % Array with all cells, [-]
fluid.ptol = 1e-4;           % Porosity tolerance to stop the simulation
fluid.Cm = 0.01;             % Injected microbial concentration, kg/m^3
fluid.Co = 0.04;             % Injected oxygen concentration, kg/m^3
fluid.Cu = 60;               % Injected urea concentration, kg/m^3

% Porosity-permeability relationship
fluid.K = @(poro) (K0.*((poro-fluid.crit)/(porosity-fluid.crit))...
        .^fluid.eta+fluid.kmin).*K0./(K0+fluid.kmin).*(poro>fluid.crit)+...
                                            fluid.kmin.*(poro<=fluid.crit);

% The two following lines are not really used in these simulations since
% the current MICP implementation only considers single-phase flow (it is
% possible to extend to two-phase flow), but since the implementation is
% based on the 'equationsOilWaterPolymer' script (two-phase flow), they are
% required to avoid errors.
fluid.bO   = fluid.bW;
fluid.rhoOS = fluid.rhoWS;

% Gravity
gravity on

%% CO2 assesment
% We simulate the CO2 distribution on the domain before MICP treatment

state0.pressure = fluid.rhoWS * norm(gravity) * c(:,3); % Initial pressure
state0.s = repmat([1, 0], G.cells.num, 1); % Initial saturations
co2     = CO2props(); % Load sampled tables of co2 fluid properties
p_ref   = 15 * mega * Pascal; % Reference pressure
t_ref   = 70 + 273.15; % Reference temperature, in Kelvin
rhoco2  = co2.rho(p_ref, t_ref); % CO2 density at ref. press/temp
cf_co2  = co2.rhoDP(p_ref, t_ref) / rhoco2; % CO2 compressibility
cf_wat  = 0; % Water compressibility (zero)
cf_rock = 0; % Rock compressibility (zero)
muco2   = co2.mu(p_ref, t_ref) * Pascal * second; % CO2 viscosity

% Use function 'initSimpleADIFluid' to make a simple fluid object
fluidH2OCO2 = initSimpleADIFluid('phases', 'WG'           , ...
                           'mu'  , [fluid.muw, muco2]     , ...
                           'rho' , [fluid.rhoWS, rhoco2]  , ...
                           'pRef', p_ref                  , ...
                           'c'   , [cf_wat, cf_co2]       , ...
                           'cR'  , cf_rock                , ...
                           'n'   , [2 2]);

% Change relperm curves
srw = 0.27;
src = 0;
fluidH2OCO2.krW = @(s) fluidH2OCO2.krW(max((s-srw)./(1-srw), 0));
fluidH2OCO2.krG = @(s) fluidH2OCO2.krG(max((s-src)./(1-src), 0));

% Add capillary pressure curve
pe = 5 * kilo * Pascal;
pcWG = @(sw) pe * sw.^(-1/2);
fluidH2OCO2.pcWG = @(sg) pcWG(max((1-sg-srw)./(1-srw), 1e-5));

% Create Model
model = TwoPhaseWaterGasModel(G, rock, fluidH2OCO2, 0, 0);

% Create Well
QCO2 = 1/day;    % Injection rate, m^3/day
r = 0.15;        % Well radius, m
cellsW = 1:1:G.cells.num;
cellsW = cellsW(c(:,1)<X(2) & c(:,3)>D+H-hb);
% Injector
WCO2 = addWell([], G, rock, cellsW, 'Type', 'rate', 'Comp_i', [0,1], ...
                                                  'Val', QCO2, 'Radius',r);

% Boundary condition
f = boundaryFaces(G);
f = f(abs(G.faces.normals(f,1)) > eps & G.faces.centroids(f,1) > X(end-1));
fp = G.faces.centroids(f,3) * fluid.rhoWS * norm(gravity);
bc = addBC([], f, 'pressure', fp, 'sat', [0 0]);

% Setup some schedule
dt = hour;
nt = 10*day/dt;
clear schedule
timestepsCO2 = repmat(dt, nt, 1);

% Make schedule
schedule = simpleSchedule(timestepsCO2, 'W', WCO2, 'bc', bc);

% Simulate
if exist('OCTAVE_VERSION', 'builtin') == 0
  fn = getPlotAfterStepCO2(state0, model, 0, 0);
  [~, statesCO2beforeMICP] = simulateScheduleAD(state0, model, schedule,...
                                                        'afterStepFn', fn);
else
  [~, statesCO2beforeMICP] = simulateScheduleAD(state0, model, schedule);
end

% Compute leakage rate (CO2 rate through the lowest grid face on the leak)
lrbeforeMICP = zeros(nt,1);
fc = G.faces.centroids;
facel =  1:1:G.faces.num;
facel = facel(fc(:,3)<D+hl+ht+.01 & fc(:,3)>D+hl+ht-.01 & ...
                                            fc(:,1)>l-a/2 & fc(:,1)<l+a/2);
for i=1:1:nt
    lrbeforeMICP(i)=abs(statesCO2beforeMICP{i}.flux(facel,2));
end

%% MICP treatment
% Create well
Q = 2.5e-4;    % Injection rate, m^3/s
Whu = 1/4;     % Upper fraction part of the well to inject the components
Whb = 1 - Whu; % Bottom fraction part of the well to inject the components
cellsW = 1:1:G.cells.num;
cellsWu = cellsW(c(:,1)<X(2) & c(:,3)>D+H-hb & c(:,3)<D+H-Whb*hb);
% Upper injector
W = addWell([], G, rock, cellsWu, 'Type', 'rate', 'Comp_i', [1,0], ...
                                                 'Val', Whu*Q, 'Radius',r);
cellsWb =cellsW(c(:,1)<X(2) & c(:,3)>D+H-Whb*hb);
% Lower injector
W = addWell(W, G, rock, cellsWb, 'Type', 'rate', 'Comp_i', [1,0], ...
                                                'Val', Whb*Q, 'Radius', r);
% The injection well is on the boundary
G.injectionwellonboundary = 1;
G.cellsinjectionwell = cellsW;

% Add the fields to the wells/bc for the additional components
bc.o = zeros(size(bc.sat,1), 1);
bc.u = zeros(size(bc.sat,1), 1);
bc.m = zeros(size(bc.sat,1), 1);
bc.b = zeros(size(bc.sat,1), 1);
bc.c = zeros(size(bc.sat,1), 1);
for i=1:2
    W(i).o = 0;
    W(i).u = 0;
    W(i).m = 0;
end
W(1).m = Cm;  % Microbes are injected since the beggining

% Create model
model = MICPModel(G, rock, fluid);

% Setup some schedule
dt = hour;
nt = 250*hour/dt;
clear schedule
timesteps = repmat(dt, nt, 1);

% Well different rates and times
N = 8; % Number of injection changes
M = zeros(N,5); % Matrix where entries per row are: time, rate, o, u, m.
M(1,1) = 15*hour/dt;
M(1,2) = Q;
M(2,1) = 26*hour/dt;
M(2,2) = eps;
M(3,1) = 100*hour/dt;
M(3,2) = Q;
M(3,3) = fluid.Co;
M(4,1) = 130*hour/dt;
M(4,2) = Q;
M(5,1) = 135*hour/dt;
M(5,2) = eps;
M(6,1) = 160*hour/dt;
M(6,2) = Q;
M(6,4) = fluid.Cu;
M(7,1) = 200*hour/dt;
M(7,2) = Q;
M(8,1) = 210*hour/dt;
M(8,2) = eps;

% Make schedule
schedule = simpleSchedule(timesteps,'W',W,'bc',bc);
for i=1:N
    schedule.control(i+1) = schedule.control(i);
    schedule.control(i+1).W(1).val = Whu*M(i,2);
    schedule.control(i+1).W(2).val = Whb*M(i,2);
    schedule.control(i+1).W(1).o = M(i,3);
    schedule.control(i+1).W(1).u = M(i,4);
    schedule.control(i+1).W(1).m = M(i,5);
    schedule.step.control(M(i,1):end) = i+1;
end

% Initial condition
state0.o = zeros(G.cells.num,1);
state0.u = zeros(G.cells.num,1);
state0.m = zeros(G.cells.num,1);
state0.b = zeros(G.cells.num,1);
state0.c = zeros(G.cells.num,1);

% If MATLAB is used, we use the getPlotAfterStepMICP function to visualize
% the results at each time step.
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    ok = 'true';
    fn = checkCloggingMICP(ok);
else
    fn = getPlotAfterStepMICP(state0, model, 0, 0);
end
[~,statesMICP] = simulateScheduleAD(state0, model, schedule, ...
                                                         'afterStepFn',fn);
%% CO2 assesment after MICP treatment
% Compute porosity and permeability after MICP treatment
porosityafterMICP = porosity - statesMICP{end}.c - statesMICP{end}.b;
KafterMICP = fluid.K(porosityafterMICP);
rock = makeRock(G, KafterMICP, porosityafterMICP);

% Create model
model = TwoPhaseWaterGasModel(G, rock, fluidH2OCO2, 0, 0);

% Make schedule
schedule = simpleSchedule(timestepsCO2, 'W', WCO2, 'bc', bc);

% Simulate
if exist('OCTAVE_VERSION', 'builtin') == 0
  fn = getPlotAfterStepCO2(state0, model, 0, 0);
  [~, statesCO2afterMICP] = simulateScheduleAD(state0, model, schedule, ...
                                                        'afterStepFn', fn);
else
  [~, statesCO2afterMICP] = simulateScheduleAD(state0, model, schedule);
end

% Compute leakage rate
lrafterMICP = zeros(max(size(timestepsCO2)),1);
for i=1:1:size(timestepsCO2)
    lrafterMICP(i)=abs(statesCO2afterMICP{i}.flux(facel,2));
end

%% Process the data
% Plot the comparison of CO2 leakage before and after micp treatment. Using
% the default example setting, we observe the leakage is reduced after MICP
% treatment but still there is significant leakage. Then more MICP
% treatments could be applied to reduce the leakage as reported in
% https://doi.org/10.1016/j.ijggc.2021.103256
figure;
hold on
plot((1:size(timestepsCO2)) * dt / day, lrbeforeMICP * 100 / QCO2, ...
                     'color', [1 .2 .2], 'LineWidth', 9, 'LineStyle', '-');
plot((1:size(timestepsCO2)) * dt / day, lrafterMICP * 100 / QCO2, ...
                      'color', [1 .5 0], 'LineWidth', 9, 'LineStyle', '-');
hold off
legend('Before MICP','After MICP','Location','best');
xlabel('Time [d]');
ylabel('CO2 leakage rate/injection rate [%]');
grid on

% If Octave is used, then the results are printed in vtk format to be
% visualize in Paraview. If MATLAB is used, then the plotToolbar is
% used to show the results.

% Write the results to be read in ParaView (GNU Octave)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mkdir vtk_2DLeakageCase;
    cd vtk_2DLeakageCase;
    mrsttovtk(G,statesCO2afterMICP, 'statesCO2afterMICP','%f');
    mrsttovtk(G,statesMICP, 'statesMICP', '%f');
    mrsttovtk(G,statesCO2beforeMICP, 'statesCO2beforeMICP','%f');
    return
end

% Plot results in GUI (MATLAB)
figure;
plotToolbar(G, statesMICP, 'field', 's:1', 'lockCaxis', true);
view(0, 0)
axis tight;
colorbar; caxis([0 1]);

%% Copyright notice
%{
Partial copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2021 NORCE Norwegian Research Centre AS, Computational
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
