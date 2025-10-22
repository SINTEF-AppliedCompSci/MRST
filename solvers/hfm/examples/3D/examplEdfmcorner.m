%% EDFM for Corner-Point Grid - Example
% This example demonstrates the use of EDFMgrid_CornerPoint for 
% non-Cartesian reservoir grids with complex fracture networks

%% Load Required Modules
mrstModule add hfm;
mrstModule add mrst-gui;
mrstModule add compositional;
mrstModule add ad-core;
mrstModule add ad-props;
mrstModule add shale;

%% Create Corner-Point Grid
nx = 20; ny = 20; nz = 10;
Lx = 200; Ly = 200; Lz = 50;

grdecl = makeModel3([nx, ny, nz], [Lx, Ly, Lz]);

% Apply elevation variations (non-Cartesian deformation)
for k = 1:nz
    factor = 1 + 0.1 * sin(2*pi*k/nz);
    startIdx = 4*(k-1)*2*nx*ny + 1;
    endIdx   = 4*k*2*nx*ny;
    grdecl.ZCORN(startIdx:endIdx) = grdecl.ZCORN(startIdx:endIdx) * factor;
end

G = processGRDECL(grdecl);
G = computeGeometry(G);

%% Define Fracture Planes
% Fracture 1: Vertical North-South
fracplanes(1).points = [50, 35, 5; 50, 180, 5; 50, 180, 45; 50, 35, 45];
fracplanes(1).aperture = 0.001;
fracplanes(1).poro = 0.5;
fracplanes(1).perm = 10000;

% Fracture 2: Vertical East-West
fracplanes(2).points = [40, 100, 10; 160, 100, 10; 160, 100, 40; 40, 100, 40];
fracplanes(2).aperture = 0.0008;
fracplanes(2).poro = 0.45;
fracplanes(2).perm = 8000;

% Fracture 3: Inclined
angle = 30 * pi/180;
fracplanes(3).points = [80, 80, 15; 120, 80, 15; ...
                        120, 120, 15 + 40*tan(angle); 80, 120, 15 + 40*tan(angle)];
fracplanes(3).aperture = 0.0005;
fracplanes(3).poro = 0.4;
fracplanes(3).perm = 5000;

% Fracture 4: Horizontal
fracplanes(4).points = [60, 60, 25; 140, 60, 25; 140, 140, 25; 60, 140, 25];
fracplanes(4).aperture = 0.0012;
fracplanes(4).poro = 0.55;
fracplanes(4).perm = 15000;
% Plot fracture planes
figure;
plotfracongrid(G,fracplanes); % visualize to check before pre-process
view(30,45)

%% Process EDFM Grid
[G, fracplanes] = EDFMgrid_CornerPoint(G, fracplanes, ...
    'Tolerance', 1e-4, ...
    'Verbose', true, ...
    'plotgrid', false);

%% Set Rock Properties
G.rock = makeRock(G, 100*milli*darcy, 0.3);

%% Compute Non-Neighboring Connections (NNCs)
tol = 1e-5;
G = fracturematrixNNC3D(G, tol);
[G, fracplanes] = fracturefractureNNCs3D(G, fracplanes, tol);

%% Define Fluid Properties
casename = 'verysimple';
[fluid, info] = getCompositionalFluidCase(casename);
eosname = 'pr';

G1cell = cartGrid([1 1], [1 1]);
G1cell = computeGeometry(G1cell);
EOSModel = EquationOfStateModel(G1cell, fluid, eosname);

p_sc = 101325;
T_sc = 288.706;
[~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);

pRef = 100*barsa;
flowfluid = initSimpleADIFluid('phases', 'WOG', 'n', [2, 2, 2], ...
                               'rho', [1000, rhoO_S, rhoG_S]);

%% Define Compositional Flow Model
gravity reset off
model = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', true);

TPFAoperators = setupShaleEDFMOpsTPFA(G, G.rock, tol);
model.operators = TPFAoperators;

%% Add Wells
totTime = 5*year;
tpv = sum(model.operators.pv);
wellRadius = 0.1;

if ~isfield(G.Matrix, 'rock')
    numMatrixCells = G.Matrix.cells.num;
    G.Matrix.rock = struct('perm', repmat(100*milli*darcy, numMatrixCells, 1), ...
                           'poro', repmat(0.3, numMatrixCells, 1));
end

centroids = G.Matrix.cells.centroids;

% Injector: Near corner (0, 0, Lz/2)
injector_point = [0, 0, Lz/2];
distances_inj = sqrt(sum((centroids - injector_point).^2, 2));
[~, sorted_idx_inj] = sort(distances_inj);
cellinj_matrix = sorted_idx_inj(1:5);

% Producer: Near corner (Lx, Ly, Lz/2)
producer_point = [Lx, Ly, Lz/2];
distances_prod = sqrt(sum((centroids - producer_point).^2, 2));
[~, sorted_idx_prod] = sort(distances_prod);
cellprod_matrix = sorted_idx_prod(1:5);

W = addWell([], G.Matrix, G.Matrix.rock, cellinj_matrix, ...
    'InnerProduct', 'ip_tpf', 'Type', 'rate', 'Val', tpv/totTime*0.7, ...
    'Radius', wellRadius, 'Comp_i', [1, 0, 0], 'Name', 'Injector', 'Sign', 1);

W = addWell(W, G.Matrix, G.Matrix.rock, cellprod_matrix, ...
    'InnerProduct', 'ip_tpf', 'Type', 'bhp', 'Val', 50*barsa, ...
    'Radius', wellRadius, 'Comp_i', [1, 1, 0]./2, 'Name', 'Producer', 'Sign', -1);

W(1).components = info.injection;
W(2).components = info.initial;

%% Visualize Grid and Wells
figure('Position', [100, 100, 1200, 600]);
subplot(1,2,1);
plotEDFMgrid(G);
hold on;
plotWell(G.Matrix, W);
view(45, 30);
title('EDFM Grid with Corner-Point Geometry');
axis tight;

%% Initialize State
s0 = [0, 1, 0];
state = initCompositionalState(G, pRef, info.temp, s0, info.initial, model.EOSModel);

%% Setup Schedule
dt = rampupTimesteps(totTime, 30*day, 8);
dt = dt(1:51);
schedule = simpleSchedule(dt, 'W', W);

%% Run Simulation
[ws, states, report] = simulateScheduleAD(state, model, schedule, ...
    'afterStepFn', getPlotAfterStep(state, model, schedule));