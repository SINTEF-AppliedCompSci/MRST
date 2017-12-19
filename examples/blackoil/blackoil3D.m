%% Water injection into a 3D fractured porous media
% Two-phase example with vertical producer/injector pair simulating
% water injection in a 3-dimensional fractured porous media using the EDFM
% method

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-core;         
mrstModule add ad-props ad-core 
mrstModule add mrst-gui;        % plotting routines
mrstModule add ad-blackoil;
mrstModule add upr;


%% NOTE ON PARALLEL PROCESSING
% Use 'parpool' to enable parallel processing. Note that on some
% machines, initializing a parallel pool takes a long time. Also note that
% MATLAB by default initializes parpool when encountering functions from
% the Parallel Processing Toolbox. To prevent this, go to Home >
% Preferences > Parallel Computing Toolbox and untick 'Automatically create
% a parallel pool...'


%% Grid and fracture(s)
% LX=500; NX=50; clx=LX/NX;
% LY=200; NY=20; cly=LY/NY;
% LZ=1; NZ=5; clz=LZ/NZ;
% celldim = [NX, NY, NZ];
% physdim = [LX, LY, LZ];
tol=1e-5;
celldim = [50, 20, 4];
physdim = [500, 200, 100];
G_matrix = cartGrid(celldim, physdim);
G_matrix = computeGeometry(G_matrix);
G_matrix.rock=makeRock(G_matrix,100*milli*darcy,0.3);

% fracplanes = struct;
% points = [30 30 0;
%           470 30 0;
%           470 30 100
%           30 30 100];
% fracplanes(1).points = points;
% fracplanes(1).aperture = 1/25;
% fracplanes(1).poro=0.8;
% fracplanes(1).perm=10000*darcy;
% fracplanes(2).points = [470 10 0
%                         470 170 0
%                         470 170 100
%                         470 10 100];
% fracplanes(2).aperture = 1/25;
% fracplanes(2).poro=0.8;
% fracplanes(2).perm=10000*darcy;
fracplanes = struct;
points = [80 100 0;
                        130 160 0;
                        130 160 100;
                        80 100 100];
fracplanes(1).points = points;
fracplanes(1).aperture = 1/25;
fracplanes(1).poro=0.8;
fracplanes(1).perm=10000*darcy;
fracplanes(2).points = [130 160 0
                        340 40 0
                        340 40 100
                        130 160 100];
fracplanes(2).aperture = 1/25;
fracplanes(2).poro=0.8;
fracplanes(2).perm=10000*darcy;
fracplanes(3).points = [250 70 0
                        330 160 0
                        330 160 100
                        250 70 100];
fracplanes(3).aperture = 1/25;
fracplanes(3).poro=0.8;
fracplanes(3).perm=10000*darcy;
% fracplanes(1).points = [250 70 0
%                         330 160 0
%                         330 160 100
%                         250 70 100];
% fracplanes(1).aperture = 1/25;
% fracplanes(1).poro=0.8;
% fracplanes(1).perm=10000*darcy;

% checkIfCoplanar(fracplanes)
plotfracongrid(G_matrix,fracplanes); % visualize to check before pre-process
G=G_matrix;


%% Process fracture(s)
[G,fracplanes]=EDFMgrid(G,fracplanes,...
    'Tolerance',tol,'plotgrid',true,'fracturelist',1:3);

%% Fracture-Matrix NNCs
G=fracturematrixNNC3D(G,tol);


%%
[G,fracplanes]=fracturefractureNNCs3D(G,fracplanes,tol);

%%
TPFAoperators = setupEDFMOperatorsTPFA(G, G.rock, tol);

%% Define fluid properties
% Define a two-phase fluid model without capillarity. The fluid model has
% the values:
%
% * densities: [rho_w, rho_o] = [1000 700 250] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 5 0.2] cP.
% * corey-coefficient: [2, 2] = [2 2 2].

fluid = initSimpleADIFluid('mu' , [   1,  5, 0.2] .* centi*poise     , ...
                           'rho', [1000, 700, 250] .* kilogram/meter^3, ...
                           'n'  , [   2,   2, 2]);

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
totTime = 5*year;
tpv = sum(model.operators.pv);
wellRadius = 0.1;
[nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
cellinj = 1:nx*ny:(1+(nz-1)*nx*ny);
cellprod1 = nx*ny:nx*ny:nz*nx*ny;
W   = addWell([], G.Matrix, G.Matrix.rock, cellinj,'InnerProduct', 'ip_tpf','Type', ...
    'rate', 'Val', tpv/totTime, 'Radius', wellRadius, 'Comp_i', [1, 0, 0], 'Name', 'Injector');
W   = addWell(W, G.Matrix, G.Matrix.rock, cellprod1, 'InnerProduct', 'ip_tpf', 'Type', ...
    'bhp' , 'Val', 50*barsa, 'Radius', wellRadius, 'Comp_i', [1, 1, 0]/2, 'Name', 'Producer1');
% cellprod2 = nx:nx*ny:(nx+(nz-1)*nx*ny);
% W   = addWell(W, G.Matrix, G.Matrix.rock, cellprod2, 'InnerProduct', 'ip_tpf', 'Type', ...
%     'bhp' , 'Val', 50*barsa, 'Radius', wellRadius, 'Comp_i', [1, 1, 0]/3, 'Name', 'Producer2');

plotWell(G,W);

%% Set up initial state and schedule
% We set up a initial state with the reference pressure and a mixture of
% water and oil initially. We also set up a simple-time step strategy that
% ramps up gradually towards 30 day time-steps.


s0 = [0.2, 0.8, 0];
state  = initResSol(G, pRef, s0);
dt = rampupTimesteps(totTime, 30*day, 8);
% dt = repmat(30*day,61,1);
schedule = simpleSchedule(dt, 'W', W);

%% Simulate problem
%fn = getPlotAfterStep(state, model, schedule);
[ws, states, report] = simulateScheduleAD(state, model, schedule, ...
   'afterStepFn', getPlotAfterStep(state, model, schedule));

%% plotting
figure;
plotToolbar(G,states);
view(40,30);
axis tight equal;

