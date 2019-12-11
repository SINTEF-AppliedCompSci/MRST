%% LOAD NECESSARY MODULES
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
mrstModule add ad-blackoil;     % AD blackoil solver
mrstModule add ad-core;         % AD core module
mrstModule add ad-props;        % AD properties

%% SET UP A STRUCTURED MATRIX GRID
physdim = [350, 200, 100]; % 350m x 200m x 100m domain
celldim = [35, 20, 10]; % 10m x 10m x 10m grid cell sizes
G = cartGrid(celldim, physdim); 
G = computeGeometry(G);
G.rock=makeRock(G,100*milli*darcy,0.3); % km=100mD, matrix porosity = 0.3

%% SET UP FRACTURE 1 (RED)
fracplanes(1).points = [40 100 0;
                        90 160 0;
                        90 160 100;
                        40 100 100]; % Vertices
fracplanes(1).aperture = 1/25;
fracplanes(1).poro=0.8;
fracplanes(1).perm=10000*darcy;

%% SET UP FRACTURE 2 (GREEN)
points = [80 160 0;
          290 40 0;
          290 40 100;
          80 160 100]; %vertices

f2normal = getnormal(points);
points([1,2],:)=points([1,2],:)-f2normal*15; % displace top points
points([3,4],:)=points([3,4],:)+f2normal*15; % displace bottom points

fracplanes(2).points = points;
fracplanes(2).aperture = 1/25;
fracplanes(2).poro=0.8;
fracplanes(2).perm=10000*darcy;

%% SET UP FRACTURE 3 (YELLOW)
fracplanes(3).points = [200 70 0;
                        280 160 0;
                        280 160 100;
                        200 70 100]; % Vertices
fracplanes(3).aperture = 1/25;
fracplanes(3).poro=0.8;
fracplanes(3).perm=10000*darcy;

%% CONSTRUCT GLOBAL GRID
[G,fracplanes]=EDFMgrid(G,fracplanes);

%% MATRIX-FRACTURE NNC CALCULATIONS
tol=1e-5; % tolerance for equality of doubles
G=fracturematrixNNC3D(G,tol);

%% FRACTURE-FRACTURE NNC CALCULATIONS
tol=1e-5; % tolerance for equality of doubles
[G,fracplanes]=fracturefractureNNCs3D(G,fracplanes,tol);

%% PLOT FRACTURE GRID
figure;
plotGrid(cartGrid([1 1 1],physdim),'facealpha',0);
hold on;
plotGrid(G,G.FracGrid.Frac1.cells.start-1+(1:G.FracGrid.Frac1.cells.num),'facecolor','r');
plotGrid(G,G.FracGrid.Frac2.cells.start-1+(1:G.FracGrid.Frac2.cells.num),'facecolor','g');
plotGrid(G,G.FracGrid.Frac3.cells.start-1+(1:G.FracGrid.Frac3.cells.num),'facecolor','y');
axis tight equal
title('Fracture Grid')
view(30,45);

%% FLUID PROPERTIES
% Define a three-phase fluid model without capillarity. Properties are
% listed in the order 'Water-Oil-Gas'. 
pRef = 100*barsa;
fluid = initSimpleADIFluid('phases' , 'WOG', ...
'mu' , [   1,  5, 0.2] * centi*poise     , ...
'rho', [1000, 700, 250] * kilogram/meter^3, ...
'c',   [1e-8, 1e-5, 1e-3] / barsa, ...
'n'  , [   2,   2, 2], ...
'pRef' , pRef);

%% DEFINE BLACK OIL MODEL
% Define a black oil model without dissolved gas or vaporized oil. Gravity
% is disabled.
gravity off
model = ThreePhaseBlackOilModel(G, G.rock, fluid, ...
'disgas', false, 'vapoil', false);
model.operators = setupEDFMOperatorsTPFA(G, G.rock, tol);

%% ADD INJECTOR
totTime = 5*year;
tpv = sum(model.operators.pv);
wellRadius = 0.1;
[nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
cellinj = 1:nx*ny:(1+(nz-1)*nx*ny);
W   = addWell([], G, G.rock, cellinj, 'Type', 'rate', ... 
'Val', tpv/totTime, 'Radius', wellRadius, ... 
'Comp_i', [1, 0, 0], 'Name', 'Injector');

%% ADD PRODUCER
cellprod = nx*ny : nx*ny : nz*nx*ny;
W   = addWell(W, G, G.rock, cellprod, 'Type', 'bhp', ... 
'Val', 50*barsa, 'Radius', wellRadius, ... 
'Comp_i', [1, 1, 0], 'Name', 'Producer');

%% INITIALIZATION
% At time zero, the model is saturated only with oil. The initial pressure
% is set to the reference pressure.
s0 = [0, 1, 0];
state  = initResSol(G, pRef, s0);

%% SET UP SCHEDULE
% Time step is set to 30 days, with an initial ramp up to the designated
% time step.
dt = rampupTimesteps(totTime, 30*day, 10);
schedule = simpleSchedule(dt, 'W', W);

%% LAUNCH SIMULATION
[ws, states, report] = simulateScheduleAD(state, model, schedule);

%% PLOT RESULTS
figure;
plotToolbar(G,states);
view(30,45);
axis equal tight
grid on