%% Simulation on the Near-wellbore modeling (NWM) hybrid grid
% 
% This example demonstrates how to generate necessary variables passed to 
% the mrst AD simulators for the NWM hybrid grid. The original data is 
% given in ECLIPSE deck file which conforms with the background 
% Corner-point grid (CPG). The class 'NearWellboreModel' accepts the data 
% for CPG and returns the variables for the hybrid grid in mrst standard 
% format, consisting of 'G', 'rock', 'fluid', 'model', 'schedule', and 
% 'initState'. Before the collection, make sure that the three subgrids 
% (GC, GV, and GW) are ready (see example 'nearWellBoreModelingGrids'). 
%
% The generation involves several key processes:
%  * Assemble the subgrdis to get the global hybrid grid
%  * Initialize the AD fluid
%  * Make rocks from subones
%  * Compute the transmissbility and neighborship (including the NNC)
%  * Setup simultaion model
%  * Convert the simulation schedule
%  * Define the initial state by equilibrium initialization
%
% Remarks:
% * In the grid domain, the subgrid boundaries are not connected. The
%   connection between subgrids are accomplished by the non-neighbor 
%   connection (NNC).
%   The 'NearWellboreModel' only accepts the ECLIPSE deck input. If you
%   need to define the simulation variables by mrst functionalities, e.g. 
%   'makeRocks', 'initSimpleADIFluid', 'addWell', and 'simpleSchedule',
%   some modifications are required.
% * This module now only support the single property and equilibration 
%   region.

mrstModule add nwm ad-core ad-blackoil ad-props mrst-gui

%% Load subgrids, well info struct and input deck 
% Load subgrids (GC, GV, GW), well info struct (well), and input deck
% (deck) in example 'nearWellBoreModelingGrids'
run nearWellBoreModelingGrids
close all

%% Define the 'NearWellboreModel'
% Define the 'NearWellboreModel' by three subgrids, input deck and well
NWM = NearWellboreModel({GC, GV, GW}, deck, well);

%% Get the global hybrid grid
G = NWM.validateGlobalGrid();

% Show the global grid. We can plot the specified subgrids by calling 
% 'G.cells.grdID==i'
figure 
subplot(1,2,1), hold on, axis off
plotGrid(G, G.cells.grdID == 1, 'facecolor', 'none')
plotGrid(G, G.cells.grdID == 2, 'facecolor', 'y')
view([-36, 38])
title('CPG and VOI grid')
subplot(1,2,2), hold on, axis off
plotGrid(G, G.cells.grdID == 2, 'facecolor', 'none')
plotGrid(G, G.cells.grdID == 3, 'facecolor', 'g')
view([-76, 60])
title('VOI grid and HW grid')
pos = get(gcf, 'position'); pos(3) = 2*pos(3);
set(gcf, 'position', pos);

% Also, use 'G.cells.grdID == i & G.cells.layers == j' to plot layer j of
% subgrid i
figure, hold on, axis off
plotCellData(G, G.cells.centroids(:,1), G.cells.layers == 4 & G.cells.grdID == 1)
plotCellData(G, G.cells.centroids(:,1), G.cells.layers == 1 & G.cells.grdID == 2)
view([-3, 56])
title('X-coordinate of the cell centroids')

%% Initialize the AD fluid
% We use 'initDeckADIFluid' to initialize the AD fluid
fluid = NWM.setupFluid();

%% Make rocks for the global grid
% First, get the rocks of three subgrids
% -------------------------------------------------------------------------
% | Rock     | Grid     | Source        | Permeability   | Anisotropy     |
% |          |          |               | coordinate     |                |
% |-----------------------------------------------------------------------|
% | rockC    | GC       | Input deck    | Local          | Yes            |
% |-----------------------------------------------------------------------|
% | rockV    | GV       | Interpolation | Global         | Yes            |
% |          |          | of rockC      |                |                |
% |-----------------------------------------------------------------------|
% | rockW    | GW       | User-defined  | Global         | No             |
% -------------------------------------------------------------------------

% Get rockC and rockV
rockC = NWM.getCPGRockFromDeck();
rockV = NWM.getVOIRocksByInterp();
% View rockC and rockV
figure
subplot(1,2,1), axis equal off
plotCellData(GC, rockC.perm(:,1), GV.parentInfo.cells{1})
title(sprintf('PermX of the CPG \nin VOI layer 1'))
subplot(1,2,2), axis equal off
plotCellData(GV, rockV.perm(:,1), GV.cells.layers==1)
title(sprintf('PermX of the VOI grid \nin VOI layer 1'))

% Define a simple rock for HW grid. Each segement (layer) of the grid
% has uniform permeability and porosity.
nclayer = GW.cells.num / GW.layers.num;
permW = linspace(400, 500, GW.layers.num) * (milli*darcy);
permW = repmat(permW, nclayer, 1);
poroW = linspace(0.18, 0.2, GW.layers.num);
poroW = repmat(poroW, nclayer, 1);
rockW.perm = [permW(:), permW(:), permW(:)];
rockW.poro = poroW(:);
% View rockW
figure, axis equal tight off, view([-85, 9])
plotCellData(GW, rockW.perm(:,1))
title('PermX of the HW grid')

% Assemble the subrocks to get the global one
rock = NWM.getGlobalRock({rockC, rockV, rockW});
% View rock
figure, hold on, axis off, view([-3, 56])
plotCellData(G, rock.perm(:,1), G.cells.layers == 4 & G.cells.grdID == 1)
plotCellData(G, rock.perm(:,1), G.cells.layers == 1 & G.cells.grdID == 2)
title('PermX of the global grid in VOI layer 1')

%% Compute the transmissbility and neighborship
% The transmissbility consists of four parts:
% -------------------------------------------------------------------------
% | Trans-     | Grid       | Flow          | Permeability   | Anisotropy |
% | missbility |            | approximation | coordinate     |            |
% |-----------------------------------------------------------------------|
% | TC         | Updated    | Linear        | Local          | Yes        |
% |            | GC         |               |                |            |
% |-----------------------------------------------------------------------|
% | TV         | Updated    | Linear        | Global         | Yes        |
% |            | GV         |               |                |            |
% |-----------------------------------------------------------------------|
% | TW         | GW         | Radial        | Global         | No         |
% |            |            |               |                |            |
% -------------------------------------------------------------------------
% | T_nnc      | Grids      | Linear        | Global         | Yes        |
% |            | connection |               |                |            |
% -------------------------------------------------------------------------
% * The T_nnc (transmissbility of NNC) is used to connect the boundaries of 
%   subgrids.
% * Updated GC and GV are obtained by 'removeCells'.

% Compute the transmissbility T of the global grid:
% T = [TC; TV; TW], corresponding to G.faces.neighbors
T = NWM.getTransGloGrid(rock);

% Generate the NNC
% First, compute the intersection relations between subgrids
intXn = NWM.computeBoundaryIntxnRelation();
% View intersection relations
figure, subplot(1,2,1), axis off
NWM.plotNonMatchingIntxnRelation(intXn ,15171)
title('Intersection relations of non-matching face: 15171')
subplot(1,2,2), axis off
NWM.plotMatchingIntxnRelation(intXn ,2689)
title('Intersection relations of matching face: 2689')
pos = get(gcf, 'position'); pos(3) = 2*pos(3);
set(gcf, 'position', pos);

% Next, generate the cell pairs and associated transmissbility
nnc = NWM.generateNonNeighborConn(intXn, rock, T);

% Get the assembled transmissbility and neighborship
% T_all = [T;                 nnc.T];
% N_all = [G.faces.neighbors; nnc.cells];
[T_all, N_all] = NWM.assembleTransNeighbors(T, nnc);

%% Setup simultaion model
% We use the 'GenericBlackOilModel' as the simulation model. The 
% 'neighbors' and 'trans' in model.operators are rewriten:
% model.operators = setupOperatorsTPFA(G, rock,'neighbors', N, 'trans', T);
% model.operators.T_all = T_all;
model = NWM.setupSimModel(rock, T_all, N_all);

%% Convert the simulation schedule
% Convert the schedule from the section 'SCHEDULE' in deck for the global
% grid. Some fields in 'W' of the HW are redefined, e.g. cells, 'WI'.
schedule = NWM.getSimSchedule(model, 'refDepthFrom', 'deck');

% Show the well
figure, axis equal tight off, view([-85, 9])
plotGrid(G, G.cells.grdID == 3, 'facecolor', 'none')
plotWell(G, schedule.control(1).W(1))

%% Get the initial state by equilibrium initialization
initState = NWM.getInitState(model);

%% Run the simulator
[wellSols, states, report] = simulateScheduleAD(initState, model, schedule);

%% Plot the results
% Oil saturation
state = states{end};
figure, axis off, view([-23, 29])
plotCellData(G, state.s(:,2))
title('Oil saturation of global grid')

figure, axis off, view([-27, 56])
plotCellData(G, state.s(:,2), G.cells.grdID==2)
title('Oil saturation of VOI grid')

figure, axis equal tight off, view([-85, 9])
plotCellData(G, state.s(:,2), G.cells.grdID==3)
title('Oil saturation of HW grid')

% Well solutions
plotWellSols(wellSols, report.ReservoirTime)