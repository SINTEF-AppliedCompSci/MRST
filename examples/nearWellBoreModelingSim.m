%% Simulation on the Near-wellbore modeling (NWM) hybrid grid
% 
% This example demonstrates how to generate necessary data structures
% passed to the mrst AD simulators for the NWM hybrid grid. The original 
% data is given in ECLIPSE deck file which conforms with the background 
% Corner-point grid (CPG). The class 'NearWellboreModel' accepts the data 
% for CPG and returns the data structures for the hybrid grid in mrst 
% standard format, consisting of 'G', 'rock', 'fluid', 'model', 'schedule', 
% and 'initState'. Before the collection, make sure that the three subgrids 
% (GC, GV, and GW) are ready (see example 'nearWellBoreModelingGrids'). 
%
% The generation involves several key processes:
%  * Assemble the subgrdis to get the global hybrid grid
%  * Initialize the AD fluid
%  * Make rocks from subones
%  * Compute the transmissibility and neighborship (including the NNC)
%  * Setup simulation model
%  * Convert the simulation schedule
%  * Define the initial state by equilibrium initialization
%
% Remarks:
% * In the grid domain, the subgrid boundaries are not connected. The
%   connections between subgrids are accomplished by the non-neighbor 
%   connection (NNC).
%   The 'NearWellboreModel' only accepts the ECLIPSE deck input. If you
%   need to define the simulation data structures by mrst functionalities,
%   e.g. 'makeRocks', 'initSimpleADIFluid', 'addWell', and 'simpleSchedule',
%   some modifications are required.
% * This module now only supports the single property and equilibration 
%   region.

clear
mrstModule add nwm ad-core ad-blackoil ad-props mrst-gui diagnostics

%% Load subgrids, well info structure and input deck 
% Load subgrids (GC, GV, GW), well info structure (well), and input deck
% (deck) in example 'nearWellBoreModelingGrids'
run nearWellBoreModelingGrids
close all

%% Define the NearWellboreModel
% Define the 'NearWellboreModel' by three subgrids, input deck and well
NWM = NearWellboreModel({GC, GV, GW}, deck, well);

%% Get the global hybrid grid
G = NWM.validateGlobalGrid();

% Show the global grid. We can plot the specified subgrids by calling 
% 'G.cells.grdID==i'
figure 
subplot(1,2,1), hold on, axis tight off
plotGrid(G, G.cells.grdID == 1, 'facecolor', 'none')
plotGrid(G, G.cells.grdID == 2, 'facecolor', 'y')
view([-36, 38])
title('CPG and VOI grid')
subplot(1,2,2), hold on, axis tight off
plotGrid(G, G.cells.grdID == 2, 'facecolor', 'none')
plotGrid(G, G.cells.grdID == 3, 'facecolor', 'g')
view([-76, 60])
title('VOI grid and HW grid')
pos = get(gcf, 'position'); pos(3) = 2*pos(3);
set(gcf, 'position', pos);

% Also, use 'G.cells.grdID == i & G.cells.layers == j' to plot layer j of
% subgrid i
figure, hold on, axis tight off
plotCellData(G, G.cells.centroids(:,1), G.cells.layers == 4 & G.cells.grdID == 1)
plotCellData(G, G.cells.centroids(:,1), G.cells.layers == 1 & G.cells.grdID == 2)
view([-3, 56])
title('X-coordinate of the cell centroids')

%% Initialize the AD fluid
% We use 'initDeckADIFluid' to initialize the AD fluid
fluid = NWM.setupFluid();

%% Make rocks for the global grid
% First, get the rocks of three subgrids:
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
subplot(1,2,1), axis equal tight off
plotCellData(GC, rockC.perm(:,1), GV.parentInfo.cells{1})
title(sprintf('PermX of the CPG \nin VOI layer 1'))
subplot(1,2,2), axis equal tight off
plotCellData(GV, rockV.perm(:,1), GV.cells.layers==1)
title(sprintf('PermX of the VOI grid \nin VOI layer 1'))

% Define a simple rock for HW grid. Each segment (layer) of the grid has 
% uniform permeability and porosity
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
figure, hold on, axis tight off, view([-3, 56])
plotCellData(G, rock.perm(:,1), G.cells.layers == 4 & G.cells.grdID == 1)
plotCellData(G, rock.perm(:,1), G.cells.layers == 1 & G.cells.grdID == 2)
title('PermX of the global grid in VOI layer 1')

%% Compute the transmissibility and neighborship
% The transmissibility consists of four parts:
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
% * The T_nnc (transmissibility of NNC) is used to connect the subgrids
% * Updated GC and GV are obtained by 'removeCells'

% Compute the transmissibility T of the global grid:
% T = [TC; TV; TW], corresponding to G.faces.neighbors
T = NWM.getTransGloGrid(rock);

% Generate the NNC
% First, compute the intersection relations between subgrids
intXn = NWM.computeIntxnRelation();
% View intersection relations
figure, subplot(1,2,1), axis off
NWM.plotNonMatchingIntxnRelation(intXn ,15171)
title('Intersection relations of non-matching face 15171')
subplot(1,2,2), axis off
NWM.plotMatchingIntxnRelation(intXn ,2689)
title('Intersection relations of matching face 2689')
pos = get(gcf, 'position'); pos(3) = 2*pos(3);
set(gcf, 'position', pos);

% Next, generate the cell pairs and associated transmissibility
nnc = NWM.generateNonNeighborConn(intXn, rock, T);

% Get the assembled transmissibility and neighborship
% T_all = [T;                 nnc.T];
% N_all = [G.faces.neighbors; nnc.cells];
[T_all, N_all] = NWM.assembleTransNeighbors(T, nnc);

%% Setup simulation model
% We use the 'GenericBlackOilModel' as the simulation model. The 
% 'neighbors' and 'trans' in model.operators are rewritten:
% model.operators = setupOperatorsTPFA(G, rock,'neighbors', N, 'trans', T);
% model.operators.T_all = T_all;
model = NWM.setupSimModel(rock, T_all, N_all);

%% Convert the simulation schedule
% Convert the schedule from the section 'SCHEDULE' in deck for the global
% grid. Some fields in 'W' of the HW are redefined, e.g. 'cells', 'WI'.
schedule = NWM.getSimSchedule(model, 'refDepthFrom', 'deck');

% Show the well
figure, axis equal tight off, view([-85, 9])
plotGrid(G, G.cells.grdID == 3, 'facecolor', 'none')
plotWell(G, schedule.control(1).W(1))

%% Get the initial state by equilibrium initialization
initState = NWM.getInitState(model);

%% Run the simulation
[wellSols, states, report] = simulateScheduleAD(initState, model, schedule);

%% Plot the results
% Oil saturation
ts = 20;
state = states{ts};

figure, axis tight off, view([-23, 29])
plotCellData(G, state.s(:,2))
title('Oil saturation of global grid')

figure, axis tight off, view([-27, 56])
plotCellData(G, state.s(:,2), G.cells.grdID==2)
title('Oil saturation of VOI subgrid')

figure, axis equal tight off, view([-85, 9])
plotCellData(G, state.s(:,2), G.cells.grdID==3)
title('Oil saturation of HW subgrid')

% Well solutions
plotWellSols(wellSols, report.ReservoirTime)

%% Apply some flow diagnostics
% Since the NNC is introduced to operators, we should rearrange the output
% flux
flux0 = state.flux;
flux0 = flux0(~all(flux0==0, 2), :);
flux = zeros(size(N_all,1), 3);
intCon = all(N_all, 2);
flux(intCon, :) = flux0;
state.flux = flux;

% Define a temporary 'G' whose 'faces.neighbors' are replaced by N_all to
% be compatible with 'computeTimeOfFlight'
Gtmp = G;
Gtmp.faces.neighbors = N_all;

% Compute TOF and tracer partitioning
W = schedule.control(2).W;
D = computeTOFandTracer(state, Gtmp, rock, 'wells', W);

% Visualize the TOF of HW subgrid
figure, axis tight off, view([-57, 67])
plotCellData(G, sum(D.tof,2), G.cells.grdID == 2)
title('TOF of HW grid')

% Visualize the injector partitioning of HW subgrid
figure, axis equal tight off, view([-85, 9])
plotCellData(G, D.ipart, G.cells.grdID == 3)
title('Injector partitioning of HW grid')