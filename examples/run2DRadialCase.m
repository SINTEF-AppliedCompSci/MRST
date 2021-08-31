%% Workflow example for ad-micp
% This example aims to show complete workflow for creating, running, and
% analyzing a 2D Radial flow system using the MICP mathematical model.

% To get distmesh for first time, uncomment and run the following 4 lines
% pth = fullfile(ROOTDIR,'utils','3rdparty','distmesh');
% mkdir(pth)
% unzip('http://persson.berkeley.edu/distmesh/distmesh.zip', pth);
% mrstPath('reregister','distmesh', pth);

% Required modules
pth = fullfile(ROOTDIR,'utils','3rdparty','distmesh');
mrstPath('reregister','distmesh', pth);
mrstModule add deckformat ad-core ad-blackoil ad-micp ad-props mrst-gui ...
                                                                   distmesh

%% Reservoir geometry/properties and model parameters
% The domain has a radius of 75 m and the grid is discretized in different
% size elements (polyhedrals). If the example is run in MATLAB, then we
% build a grid using the cylinder function.

% Grid
R = 75;       % Reservoir radius, m

if exist('OCTAVE_VERSION', 'builtin') ~= 0 % GNU Octave
    % The grid is coarser in order to run the example in couple of minutes.
    % You should try with finer grids (e.g., hmin=.75 and hmax=10).
    B = 20;       % hmin to hmax transition radius, m
    hmin = 4;     % Minimum grid size, m
    hmax = 30;    % Maximum grid size, m

    fd = @(p) dcircle(p,0,0,R);
    fh = @(p) hmin+.3*abs(dcircle(p,0,0,0)).* (abs(dcircle(p,0,0,0))<B)+...
      min(hmin+.3*abs(dcircle(p,0,0,B)),hmax).* (abs(dcircle(p,0,0,0))>=B);
    [p,t] = distmesh2d(fd, fh, hmin, [-R,-R;R,R], ...
                                                [-R,-R;R,-R;-R,R;R,R;0,0]);
    G = makeLayeredGrid(pebi(triangleGrid(p, t)),1);
    close
else % MATLAB
    nc = 20;                  % Number of partitions around the circle
    P = [];
    rr = -1.25:0.25:0;
    hmin = 1;                                   % Minimum grid size, m
    hmax = R*(exp(rr(end))-exp(rr(end-1)));     % Maximum grid size, m
    for r = [hmin:3:18 R*exp(rr)]
        [x,y,z] = cylinder(r,nc);
    P = [P [x(1,:); y(1,:)]];
    end
    P = unique([P'; 0 0],'rows');
    G = makeLayeredGrid(pebi(triangleGrid(P)),1);
end
G = computeGeometry(G);
c = G.cells.centroids;
C = ones(G.cells.num,1);

% Rock
K0 = 1e-12*C;               % Aquifer permeability, m^2
porosity = 0.2;             % Aquifer porosity, [-]
rock = makeRock(G, K0, porosity);

% Fluid properties
fluid.muw = 2.535e-4;        % Water viscocity, Pa s
fluid.bW   =  @(p) 0*p + 1;  % Water formation volume factor, [-]
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
%% Define well, boundary conditions, and simulation schedule
% The components are injected from the middle point of the reservoir and
% the free-flow boundary on the boundary domain is approximated by constant
% pressure.

% Create Well
Q = 1.2e-3;   % Injection rate, m^3/s
r = 0.15;     % Well radius, m
cellsW =  G.cells.indexMap;
cellsW = cellsW(c(:,1).^2 + c(:,2).^2<(hmin/2)^2);
% Injector
W = addWell([], G, rock, cellsW, 'Type', 'rate', 'Comp_i', [1,0], ...
                                                    'Val', Q, 'Radius', r);
% Add the fields to the wells for the additional injected components
W.o = 0;
W.u = 0;
W.m = fluid.Cm;   % Microbes are injected since the beggining

% The injection well is not on the boundary
G.injectionwellonboundary = 0;

% Boundary condition
f = boundaryFaces(G);
f = f((abs(G.faces.normals(f,2))>eps | abs(G.faces.normals(f,1))>eps) ...
           & G.faces.centroids(f,1).^2+G.faces.centroids(f,2).^2>R-hmax/4);
bc = addBC([], f, 'pressure', atm, 'sat', [0 0]);
bc.o = zeros(size(bc.sat,1), 1);
bc.u = zeros(size(bc.sat,1), 1);
bc.m = zeros(size(bc.sat,1), 1);
bc.b = zeros(size(bc.sat,1), 1);
bc.c = zeros(size(bc.sat,1), 1);

% Setup some schedule
dt = hour;
nt = 300*hour/dt;
clear schedule
timesteps = repmat(dt, nt, 1);

% Well different tates and times
N = 8; % Number of injection changes
M = zeros(N,5); % Matrix where entries per row are:time, rate, o, u, m.
M(1,1) = 20*hour/dt;
M(1,2) = Q;
M(2,1) = 40*hour/dt;
M(2,2) = eps;
M(3,1) = 140*hour/dt;
M(3,2) = Q;
M(3,3) = fluid.Co;
M(4,1) = 160*hour/dt;
M(4,2) = Q;
M(5,1) = 180*hour/dt;
M(5,2) = eps;
M(6,1) = 230*hour/dt;
M(6,2) = Q;
M(6,4) = fluid.Cu;
M(7,1) = 250*hour/dt;
M(7,2) = Q;
M(8,1) = 270*hour/dt;
M(8,2) = eps;

% Make schedule
schedule = simpleSchedule(timesteps,'W',W,'bc',bc);
for i=1:N
    schedule.control(i+1) = schedule.control(i);
    schedule.control(i+1).W(1).val = M(i,2);
    schedule.control(i+1).W(1).o = M(i,3);
    schedule.control(i+1).W(1).u = M(i,4);
    schedule.control(i+1).W(1).m = M(i,5);
    schedule.step.control(M(i,1):end) = i+1;
end

%% Set up simulation model
model = MICPModel(G, rock, fluid);

%% Define initial state
% Initially, there is only water without dissolved components.

state0      = initState(G, W, atm, [1, 0]);
state0.o    = zeros(G.cells.num,1);
state0.u    = zeros(G.cells.num,1);
state0.m    = zeros(G.cells.num,1);
state0.b    = zeros(G.cells.num,1);
state0.c    = zeros(G.cells.num,1);

%% Simulate 2DRadialCase
% If MATLAB is used, we use the getPlotAfterStepMICP function to visualize
% the results at each time step.

% Simulate case (GNU Octave/MATLAB)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    ok = 'true';
    fn = checkCloggingMICP(ok);
else
    fn = getPlotAfterStepMICP(state0, model, 0, 90);
end
[~, states] = simulateScheduleAD(state0, model, schedule,'afterStepFn',fn);

%% Process the data
% If Octave is used, then the results are printed in vtk format to be
% visualize in Paraview. If MATLAB is used, then the plotToolbar is
% used to show the results.

% Write the results to be read in ParaView (GNU Octave)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mkdir vtk_2DRadialCase;
    cd vtk_2DRadialCase;
    mrsttovtk(G,states,'states','%f');
    return
end

% Plot results in GUI (MATLAB)
figure;
plotToolbar(G, states,'field', 's:1','lockCaxis',true);
view(0,90)
axis tight;
colorbar; caxis([0 1]);

%% Copyright notice
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
