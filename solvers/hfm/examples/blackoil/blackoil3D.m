%% Water injection into a 3D fractured porous media
% Two-phase example with vertical producer/injector pair simulating
% water injection in a 3-dimensional fractured porous media using the EDFM
% method

%% Load necessary modules
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
mrstModule add ad-blackoil;     % AD blackoil solver
mrstModule add ad-core;         % AD core module
mrstModule add ad-props;        % AD properties

%% Set up a matrix grid
% We first set up a Cartesian matrix grid with dimensions 500m x 200m x
% 100m with grid blocks that are 10m x 10m x 10m. Matrix permeability is
% 100mD and porosity is 30%.

celldim = [50, 20, 10];
physdim = [500, 200, 100];
G = cartGrid(celldim, physdim); 
G = computeGeometry(G);
G.rock=makeRock(G,100*milli*darcy,0.3);

% Plot Matrix Grid
figure;
plotGrid(G,'facealpha',0,'edgealpha',0.5);
view(30,45);

%% Set up fracture planes
% Fracture planes are set up by defining their vertices. Additionally,
% the aperture, porosity and permeability of the fractures are provided.
% Fracture planes 1 and 3 will be vertical while fracture 2 is slanted.

% Fracture plane 1
fracplanes(1).points = [90 100 0;
                        140 160 0;
                        140 160 100;
                        90 100 100];
fracplanes(1).aperture = 1/25;
fracplanes(1).poro=0.8;
fracplanes(1).perm=10000*darcy;

% Fracture plane 2 
points = [130 160 0
          340 40 0
          340 40 100
          130 160 100];

f2normal = getnormal(points);
points([1,2],:)=points([1,2],:)+f2normal*5; % displace top points
points([3,4],:)=points([3,4],:)-f2normal*5; % displace bottom points

fracplanes(2).points = points;
fracplanes(2).aperture = 1/25;
fracplanes(2).poro=0.8;
fracplanes(2).perm=10000*darcy;

% Fracture plane 3
fracplanes(3).points = [250 70 0
                        330 160 0
                        330 160 100
                        250 70 100];
fracplanes(3).aperture = 1/25;
fracplanes(3).poro=0.8;
fracplanes(3).perm=10000*darcy;

% Plot fracture planes
figure;
plotfracongrid(G,fracplanes); % visualize to check before pre-process
view(30,45)

%% Construct fracture grid
% The fracture grid is constructed using the matrix grid. The matrix grid
% will serve as a 'cookie cutter' to subdivide the fracture planes.
% Parallel processing can be used to speed up this process (Start a 
% parallel pool to do this).

tol=1e-5;
[G,fracplanes]=EDFMgrid(G,fracplanes,...
    'Tolerance',tol,'fracturelist',1:3);

% Plot Fracture grid
figure;
plotGrid(cartGrid([1 1 1],physdim),'facealpha',0);
hold on;
plotGrid(G,G.Matrix.cells.num+1:G.cells.num);
axis tight equal
title('Fracture Grid')
view(30,45);

%% Fracture-Matrix non-Neighbouring Connections (NNC)
% This calculates the transmissibilities between connected matrix and
% fracture grid blocks. Information is saved under G.nnc.

tol=1e-5;
G=fracturematrixNNC3D(G,tol);

% Plot matrix gridblocks that have NNCs with fractures
figure;
plotfracongrid(cartGrid([1 1 1],physdim),fracplanes);
hold on;
plotGrid(G,G.nnc.cells(:,1),'facealpha',0,'edgealpha',1);
axis tight equal
title('Matrix-Fracture NNCs')
view(30,45);

%% Fracture-Fracture NNCs
% This calculates the transmissibilities between connected fracture and
% fracture grid blocks. Information is saved under G.nnc.

tol=1e-5;
[G,fracplanes]=fracturefractureNNCs3D(G,fracplanes,tol,'Verbose',true);

% Plot fracture gridblocks that have Fracture-Fracture NNCs
fraccells = G.nnc.cells(strcmp(G.nnc.type,'fracfrac'),:);
figure;
plotGrid(cartGrid([1 1 1],physdim),'facealpha',0);
hold on;
plotGrid(G,G.Matrix.cells.num+1:G.cells.num,'facealpha',0,'edgealpha',0.5);
plotGrid(G,fraccells);
axis equal tight;
title('Fracture-Fracture NNCs')
view(30,45);


%% Setup TPFA Operators
% Generate operators for the black oil model which incorporates all the
% NNCs identified above.

TPFAoperators = setupEDFMOperatorsTPFA(G, G.rock, tol);

%% Define fluid properties
% Define a three-phase fluid model without capillarity.

fluid = initSimpleADIFluid('mu' , [   1,  5, 0.2] .* centi*poise     , ...
                           'rho', [1000, 700, 250] .* kilogram/meter^3, ...
                           'n'  , [   2,   2, 2]);

% Add compressibility to fluid
pRef = 100*barsa;
c_w = 1e-8/barsa;
c_o = 1e-5/barsa;
c_g = 1e-3/barsa;

fluid.bW = @(p) exp((p - pRef)*c_w);
fluid.bO = @(p) exp((p - pRef)*c_o);
fluid.bG = @(p) exp((p - pRef)*c_g);

%% Define three-phase compressible flow model
% We define a three-phase black-oil model without dissolved gas or vaporized
% oil. This is done by first instantiating the blackoil model, and then
% manually passing in the internal transmissibilities and the topological
% neighborship from the embedded fracture grid.

gravity reset off
model = ThreePhaseBlackOilModel(G, [], fluid, 'disgas', false, 'vapoil', false);
model.operators = TPFAoperators;

%% Add wells
% An injector/producer pair is added. 1PV of water is injected over the
% course of 5 years.

totTime = 5*year;
tpv = sum(model.operators.pv);
wellRadius = 0.1;
[nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));

% Injector
cellinj = 1:nx*ny:(1+(nz-1)*nx*ny);
W = addWell([], G.Matrix, G.Matrix.rock, cellinj, ...
            'InnerProduct', 'ip_tpf','Type',  'rate', 'Val', ...
            tpv/totTime, 'Radius', wellRadius, 'Comp_i', [1, 0, 0], ...
            'Name', 'Injector', 'Sign', 1);

% Producer
cellprod = nx*ny:nx*ny:nz*nx*ny;
W = addWell(W, G.Matrix, G.Matrix.rock, cellprod, ...
            'InnerProduct', 'ip_tpf', 'Type', 'bhp', 'Val', 50*barsa, ...
            'Radius', wellRadius, 'Comp_i', [1, 1, 0], ...
            'Name', 'Producer', 'Sign', -1);

% Plot setup
plotEDFMgrid(G);
hold on;
plotWell(G,W);

%% Set up initial state and schedule
% We set up a initial state with the reference pressure and a mixture of
% water and oil initially. We also set up a simple-time step strategy that
% ramps up gradually towards 30 day time-steps.

s0 = [0.2, 0.8, 0];
state  = initResSol(G, pRef, s0);
dt = rampupTimesteps(totTime, 30*day, 8);
schedule = simpleSchedule(dt, 'W', W);

%% Simulate problem
[ws, states, report] = simulateScheduleAD(state, model, schedule, ...
   'afterStepFn', getPlotAfterStep(state, model, schedule));
