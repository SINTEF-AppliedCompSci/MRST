clearvars -except METISPATH mrstVerbose screenSize
min_size = 1; cell_size = 5;
a = 1/25; % fracture aperture
dof_frac = 1; % Coarse dof
dof_matrix = 50; % Coarse dof
layers = 5;
%% User defined fracture lines and appropriate grid

nx = 100; ny = 100;
G = cartGrid([nx ny]);
G = computeGeometry(G);
fl = [30 30 70 70; 30 70 70 30];
% load('Lines13_Network1_100by100Domain.mat'); nx = G.cartDims(1); ny = G.cartDims(2);
flayers = 2:4;


%% Process fracture lines into grid

[G,fracture] = processFracture2D(G,fl); fracture.aperture = a;
% figure; plotFractureLines(G,fracture,'network'); % 3rd arg can be line #'s, network or none
% figure; plotMarkedCells(G, fracture,'lines'); % 3rd arg can be line #'s, network or none


%% Compute CI and grid fracture

dispif(mrstVerbose, 'Computing CI and gridding fracture...\n\n');
G = CIcalculator2D(G,fracture);
% plotCI(G,fracture);

[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
% figure; plotFractureNodes2D(G,F,fracture,show_node_num); clear show_node_num


%% Make Layered Grid
dispif(mrstVerbose, 'Extruding...\n\n');
Gl = makeLayers(G,flayers,layers);


%% Rock and fluid properties
Gl.rock.perm = ones(Gl.cells.num,1)*milli*darcy();
Gl.rock.poro = 0.5*ones(Gl.cells.num, 1);
K_frac = 1000; % Scaling factor = K_frac/K_mat D
poro_frac = 0.5;
Gl = makeRockFrac(Gl, K_frac, 'permtype','homogeneous','rockporo',poro_frac);

%% Global grid with NNC's and corresponding transmissibility
dispif(mrstVerbose, 'Define NNC''s...\n\n');
[Gl,T] = makeNNCgrid(G,Gl,F,fracture,flayers);
G = Gl;


%% Init

pinit = 300*barsa();
state = initResSol (G, pinit, 0);

radius = 1e-2;
cellinj = 1:nx*ny:nx*ny*layers;
cellprod = nx*ny:nx*ny:nx*ny*layers;
W   = addWell([], G.Matrix, G.Matrix.rock, cellinj, 'Type', 'bhp',...
    'Val', 350*barsa, 'Radius', radius, 'Sign',1, 'Comp_i', [1, 0]);
W   = addWell(W, G.Matrix, G.Matrix.rock, cellprod, ...
    'Type', 'bhp', 'Val', 250*barsa(), 'Radius', radius, 'Sign', -1, 'Comp_i', [0, 1]);

bc = [];
state = incompTPFA(state, G, T, fluid,  ...
    'bc',bc, 'Wells', W, 'MatrixOutput', true, 'use_trans',true);

nt     = 30;
Time   = 100*year;
dT     = Time/nt;


%% Set up AD model
mrstModule add ad-core ad-blackoil blackoil-sequential ad-props

fluidAD = initSimpleADIFluid('rho', [1000, 800, 0], 'mu', [1, 1, 1]*centi*poise, 'n', [1 1 1]);

% Add some slight compressibility
c = 0.001/barsa;
p_ref = 300*barsa;
fluidAD.bO = @(p) exp((p - p_ref)*c);

% Send in empty rock to avoid setting up operators as it does not account
% for NNC or alternate trans
model = TwoPhaseOilWaterModel(G, [], fluidAD);
N = getNeighbourship(G, 'topological', true);
intx = all(N ~= 0, 2);
% Send in internal transmissibility and neighborship
model.operators = setupOperatorsTPFA(G, G.rock, 'trans', T(intx), 'neighbors', N(intx, :));

%% Set up schedule
timesteps = repmat(dT, nt, 1);
schedule = simpleSchedule(timesteps, 'Wells', W);

state0 = initResSol (G, pinit, [0, 1]);

%% Fully implicit with preconditioner
linearsolver = CPRSolverAD();
[wsFI, statesFI, repFI] = simulateScheduleAD(state0, model, schedule, 'linearsolver', linearsolver);

%% Sequential pressure/transport
modelSeq = getSequentialModelFromFI(model);
[ws, states, rep] = simulateScheduleAD(state0, modelSeq, schedule);


%% Dynamic timestepping
scheduleBig = compressSchedule(schedule);

solver = NonLinearSolver();
solver.timeStepSelector = IterationCountTimeStepSelector('firstRampupStep', 30*day);
solver.timeStepSelector.targetIterationCount = 3;

[wsAdaptive, statesAdaptive, repAdaptive] = simulateScheduleAD(state0, modelSeq, scheduleBig, 'NonLinearSolver', solver, 'outputministeps', true);

[~, timesteps] = convertReportToSchedule(repAdaptive, scheduleBig);

%% Compare well sols
Tstep =  {cumsum(schedule.step.val), cumsum(timesteps), cumsum(schedule.step.val)};
wellsols = {wsFI, wsAdaptive, ws};
names = {'Fully implicit', 'Sequential implicit (adaptive)', 'Sequential implicit'};

plotWellSols(wellsols, Tstep, 'datasetnames', names);