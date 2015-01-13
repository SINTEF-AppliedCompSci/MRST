%% POLYMER BOUNDARY CONDITIONS AND SOURCE EXAMPLE
%
% This example is made just to illustrate how one can setup a polymer
% simulation with boundary conditions and/or source.
% 

% Required modules
mrstModule add deckformat ad-core ad-fi mrst-gui


%% Setup case

% Grid, rock and fluid
deck  = readEclipseDeck('POLYMER.DATA');
deck  = convertDeckUnits(deck);
G     = initEclipseGrid(deck);
G     = computeGeometry(G);
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
fluid = initDeckADIFluid(deck);

% Gravity
gravity on

% Initial state
state0      = initResSol(G, 100*barsa, [ .1, .9]);
state0.c    = zeros(G.cells.num,1);
state0.cmax = zeros(G.cells.num,1);

% Create model
model = OilWaterPolymerModel(G, rock, fluid);

% Setup some schedule
dt = 25*day;
nt = 10;
clear schedule
schedule.step.val       = dt.*ones(nt,1);
schedule.step.control   = ones(nt,1);
schedule.control(1).W   = [];
schedule.control(1).bc  = [];
schedule.control(1).src = [];


%% Pressure (Dirichlet) Boundary Condition

% Create Dirichlet boundary condition
bc = pside([], G, 'xmin', 300*barsa, 'sat', [1 0]);
bc = pside(bc, G, 'xmax', 100*barsa, 'sat', [0 0]);
bc.poly = 4.*ones(size(bc.sat,1), 1);
schedule.control(1).bc  = bc;
schedule.control(1).src = [];

% Simulate
[~, states] = simulateScheduleAD(state0, model, schedule);

% Plot results in GUI
figure;
plotToolbar(G, states);
view([-10, 14]);
axis tight;
colorbar;

%% Flux (Neumann) Boundary Condition

% Create Neumann boundary condition
bc = fluxside([], G, 'xmin',  0.004, 'sat', [1 0]);
bc = fluxside(bc, G, 'xmax', -0.004, 'sat', [0 0]);
bc.poly = 4.*ones(size(bc.sat,1), 1);
schedule.control(1).bc  = bc;
schedule.control(1).src = [];

% Simulate
[~, states] = simulateScheduleAD(state0, model, schedule);

% Plot results in GUI
figure;
plotToolbar(G, states);
view([-10, 14]);
axis tight;
colorbar;


%% Source

% Create source
ijk = gridLogicalIndices(G);
srcCells = find(ijk{1}==5  & ijk{2}==5  & ijk{3}==2);
snkCells = find(ijk{1}==26 & ijk{2}==26 & ijk{3}==2);
srcVals  = 0.004.*ones(numel(srcCells),1);
src = addSource( [], srcCells,  srcVals, 'sat', [1 0]);
src = addSource(src, snkCells, -srcVals, 'sat', [1 0]);
src.poly = 4.*ones(size(src.sat,1), 1);
schedule.control(1).bc  = [];
schedule.control(1).src = src;

% Simulate
[~, states] = simulateScheduleAD(state0, model, schedule);

% Plot results in GUI
figure;
plotToolbar(G, states);
view([-10, 14]);
axis tight;
colorbar;

