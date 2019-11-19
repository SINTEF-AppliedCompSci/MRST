%% Near-wellbore model (NWM) coupling with the Multi-Segment well (MSW) model 
% This model shows how to couple the NWM with MSW. The 'MultiSegWellNWM'
% will automatically generate the nodes and segments for the MSW simulation
% through building a 'wellbore grid'. The nodes of MSW correspond to the
% cells of wellbore grid, and the segments of MSW correspond to the faces
% of wellbore grid. The geometrical information can be obtained by grid
% geometries, e.g. 'cells.centroids', 'cells.volumes'. The topology of the
% segments is 'faces.neighbors'.

clear
% Load necessary modules
mrstModule add ad-core ad-blackoil ad-props mrst-gui wellpaths deckformat
%% Read the ECLIPSE input deck
fn = fullfile(pwd, '..', 'data', 'CPG.data');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);
%% Build the background Corner-point grid (CPG)
GC = initEclipseGrid(deck);
GC = computeGeometry(GC);
%% Define basic information of horizontal well (HW)
% Define the well trajectory
load trajectory.mat

% Number of well segments
ns = size(pW,1)-1;

% Define the well struct. 
%  'name':         Name of the well, should match the well name list in deck 
%  'trajectory':   Well trajectory, 3D points in xyz format
%  'segmentNum':   Number of well segment, equals to n_wellpoints - 1
%  'radius':       Wellbore (Casing) radius per well segment
%  'skinFactor':   Skin factor per well node
%  'openedSegs':   Opened segments (allow fluid to flow into wellbore)
% The coupling of the NWM with multi-segment well requires: 
%  ** 'isMS'      :   Indicating the multi-segment well definition
%  ** 'roughness' :   Roughness per well segment

well = struct(...
    'name'         , 'P1', ...
    'trajectory'   , pW, ...
    'segmentNum'   , ns, ...
    'radius'       , 0.15 * ones(ns+1,1), ...
    'skinFactor'   , zeros(ns,1), ...
    'openedSegs'   , 1:ns, ...
    'isMS'         , true, ...
    'roughness'    , 0 * ones(ns,1));
%% Define the volume of interest (VOI)
% Define the 2D boundary of VOI. The well should be located inside this 
% boundary in xy plane.
pbdy = [240,   50;...
        160,   80;...
        120,  160;...
        150,  205;...
        230,  170;...
        280,   90];

% Define the extra layers:
nextra = [1, 1];

% Define the VOI according to the CPG, well, boundary and extra layers
VOI = VolumeOfInterest(GC, well, pbdy, nextra);

% Get the geometrical information of CPG in VOI, including cells, faces,
% boundary faces, nodes, boundary nodes, etc.
packed = allInfoOfVolume(VOI);
%%
close all
%% Build the layered unstructured VOI grid
% The number of Cartesian cells in Y direction
WR.ny = 6; 
% The size Cartesian region in Y direction, better to be larger than
% the well-segment length
VOI.maxWellSegLength2D()
WR.ly = 12;
% The number of angular cells in radial region
WR.na = 5;
% Prepare the 2D well region nodes
WR = VOI.prepareWellRegionNodes2D(WR);

% Plot the 2D well region grid. This grid will be glued to the unstructured
% grid.
VOI.plot2DWRSubGrid(WR)

% Define the number of refined layers for each VOI layer. Each of the four 
% VOI layers will be refined into two layers.
layerRf = [2, 2, 2, 2];

% Reconstruct the CPG in VOI to layered unstructured grid
GV = VOI.ReConstructToUnstructuredGrid(WR, layerRf, ...
    'multiplier', 0.2, 'maxIter', 500, 'gridType', 'Voronoi');

% Show the VOI grid
figure, axis off
arrayfun(@(layer)plotGrid(GV, GV.cells.layers==layer, ...
    'facecolor', rand(3,1)), 1:GV.layers.num)
view([-41, 62])
%%
close all
%% Build the layered radial HW grid

%               [ymin, ymax, zmin, zmax]
regionIndices = [   2,    5,    2,    7];

% Define the HW region according to the GV, well, and regionIndices
HW = HorWellRegion(GV, well, regionIndices);

% Visualize the HW region
HW.showWellRegionInVOIGrid('showWellRgionCells', true);
view([-53, 15]), axis off

% Get the geometrical information of GV in HW region, including cells, 
% faces, boundary faces, nodes, boundary nodes, etc.
packedW = HW.allInfoOfRegion();

radPara1 = struct(...
    'gridType'  , 'pureCircular', ...
    'maxRadius' , 2, ...
    'nRadCells' , 8);

radPara2 = struct(...
    'gridType'  , 'gradual', ...
    'boxRatio'  , [0.7, 0.7], ...
    'nRadCells' , [7, 2], ...
    'pDMult'    , 10, ...
    'offCenter' , true);

% Reconstruct the VOI grid in HW region to layered radial grid
GW = HW.ReConstructToRadialGrid(radPara2);

% Visualize the raidal grid in HW region
HW.showWellRegionInVOIGrid('showWellRgionCells', false);
plotGrid(GW, GW.cells.layers==1, 'facecolor', rand(3,1))
view([-53, 15]), axis off

% Show the HW grid
figure, axis equal tight off, view([-85, 9])
arrayfun(@(layer)plotGrid(GW, GW.cells.layers==layer, ...
    'facecolor', rand(3,1)), 1:GW.layers.num)
%%
close all
%% Collect the simulation data
% Define the Multi-Segment well NearWellboreModel by three subgrids, input 
% deck and well
MSW = MultiSegWellNWM({GC, GV, GW}, deck, well);
%% Get the global hybrid grid
G = MSW.validateGlobalGrid();
%%
close all
%% Setup the fuild
fluid = MSW.setupFluid();
%% Get the rocks
% CPG rock
rockC = MSW.getCPGRockFromDeck();
% VOI grid rock
rockV = MSW.getVOIRocksByInterp();

% Define a simple rock for HW grid. Each segement (layer) of the grid
% has uniform permeability and porosity.
nclayer = GW.cells.num / GW.layers.num;
permW = linspace(400, 500, GW.layers.num) * (milli*darcy);
permW = repmat(permW, nclayer, 1);
poroW = linspace(0.18, 0.2, GW.layers.num);
poroW = repmat(poroW, nclayer, 1);
rockW.perm = [permW(:), permW(:), permW(:)];
rockW.poro = poroW(:);

% Get the rock for global grid
rock = MSW.getGlobalRock({rockC, rockV, rockW});
%%
close all
%% Get the transmissbility and neighbors
% Compute the transmissbility of the global grid
T = MSW.getTransGloGrid(rock);

% Generate NNC
intXn = MSW.computeBoundaryIntxnRelation();
nnc = MSW.generateNonNeighborConn(intXn, rock, T);

% Get the assembled transmissbility and neighbors
[T_all, N_all] = MSW.assembleTransNeighbors(T, nnc);
%%
close all
%% Setup simultaion model
model = MSW.setupSimModel(rock, T_all, N_all);

%% Get the simulation schedule
% The nodes and segments for the MSW simulation are generated through 
% building a 'wellbore grid'. The nodes of MSW correspond to the cells of
% wellbore grid, and the segments of MSW correspond to the faces of 
% wellbore grid. The geometrical information can be obtained by grid 
% geometries, e.g. 'cells.centroids', 'cells.volumes'. The topology of the
% segments is 'faces.neighbors'.
gW  = MSW.wellboreGrid;
nodes = MSW.generateNodes()
segs = MSW.generateSegments()

% Plot the reservoir cells connected to node 2
MSW.plotCell2Node(nodes, 2)
axis equal tight off, view([-85, 9])

% Define the schedule (the reference depth has been reset to the depth of
% top node)
schedule = MSW.getSimSchedule(model);

% Show the well
figure, axis equal tight off, view([-85, 9])
plotGrid(G, G.cells.grdID == 3, 'facecolor', 'none')
plotWell(G, schedule.control(1).W(3))
%% Get the initial state by equilibrium initialization
initState = MSW.getInitState(model);
%% Simulate base case
[wellSols, states, report] = simulateScheduleAD(initState, model, schedule);
%% Plot the results
state = states{end};

% Oil saturation
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
%% Plot the node pressure 
ws = wellSols{end}(3);
nPres = [ws.bhp; ws.nodePressure];

figure, axis equal tight off, view([-85, 9])
plotCellData(G, state.pressure, G.cells.grdID==3 & G.cells.layers==1)
plotCellData(gW, nPres, gW.cells.layers==1)

W = schedule.control(1).W(3);
L = cumsum([0; W.segments.length]);
figure
plot(L, nPres/barsa, 's-')
xlabel('Distance from top node (m)')
ylabel('Node pressure (barsa)')
