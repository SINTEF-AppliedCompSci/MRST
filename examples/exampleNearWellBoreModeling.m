%% Near-wellbore modeling (NWM) example
% 
% Typically, the horizontal well (HW) is built in the Corner-point grid 
% (CPG) or Cartesian grid by the standard well model (Peaceman's model) 
% that inserts the virtual well into coarse reservoir grid. However, the 
% standard model gives a grid-based well trajectory, which will be rather 
% unrigorous when the real well trajectory fails to follow the grid 
% orientation. For another, some field operations require high-resolution 
% flow descriptions in the well vicinity. The coarse well cells always 
% unable to provide such descriptions due to the linear flow approxmation 
% and low-resolution rock properties. To this end, this example 
% demonstrates the local grid reconstruction in the near-wellbore region of
% horizontal well. 
%
% The global grid will consist of three subgrids that adopt different
% gridding strategies and transmissibility calculations:
% -------------------------------------------------------------------------
% | Grid | Description | Type                 | Transmissibility          |
% |-----------------------------------------------------------------------|
% | GC   | Background  | Corner-point         | Local coordinate system   |
% |      | grid        |                      | Linear approximation      |
% |-----------------------------------------------------------------------|
% | GV   | VOI grid    | Unstructured,        | Global coordinate system  |
% |      |             | Vertically layered   | Linear approximation      |
% |-----------------------------------------------------------------------|
% | GW   | HW grid     | Structured radial,   | Global coordinate system  |
% |      |             | Horizontally layered | Radial approximation      |
% |      |             |                      | Isotropic permeability    |
% -------------------------------------------------------------------------
%
% Remarks:
% * In the grid domain, the subgrid boundaries are not connected. The
%   connection between subgrids are accomplished by the non-neighbor 
%   connection (NNC).
% * The ECLIPSE input deck is used to specify basic simulation parameters,
%   including 'RUNSPEC', 'GRID', 'PROPS', 'SOLUTION', and 'SCHEDULE'. 
%   These parameters can also be provided by MRST functionalities, e.g. 
%   'makemodel3', 'processGRDECL', 'initSimpleADIFluid', 'addWell', 
%   'simpleSchedule'. However, this module is basically tailored to the 
%   deck-format input. Some modifications of this module's function are 
%   needed when calling above functions.
% * This module now only support the modeling of horizontal wells.
% * This module now only support the single property and equilibration 
%   region.
% * The volume of interest (VOI) can not cover the faults.
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
% The well trajectory is specified by a set of discrete 3D well points 
% (in xyz format) which divides the HW into multiple segments.
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
%  'isMS'      :   Indicating the multi-segment well definition
%  'roughness' :   Roughness per well segment
% (See 'exampleNWMMultiSegmentWell')

well = struct(...
    'name'         , 'P1', ...
    'trajectory'   , pW, ...
    'segmentNum'   , ns, ...
    'radius'       , 0.15 * ones(ns+1,1), ...
    'skinFactor'   , zeros(ns,1), ...
    'openedSegs'   , (1:ns));
%% Define the volume of interest (VOI)
% The VOI is a 3D region in which we will reconstruct a layered 
% unstructured grid. 
% Define the 2D boundary of VOI. The well should be located inside this 
% boundary in xy plane.
pbdy = [240,   50;...
        160,   80;...
        120,  160;...
        150,  205;...
        230,  170;...
        280,   90];

% The VOI are vertically expanded by extra layers:
% nextra(1): Layer numbers above the layers that HW occupies
% nextra(2): Layer numbers below the layers that HW occupies
nextra = [1, 1];

% Define the VOI according to the CPG, well, boundary and extra layers
VOI = VolumeOfInterest(GC, well, pbdy, nextra);

% Get the geometrical information of CPG in VOI, including cells, faces,
% boundary faces, nodes, boundary nodes, etc.
geoV = VOI.allInfoOfVolume();

% Show the VOI boundary, cells, and faces.
VOI.plotVolumeCells(geoV), view(3)
VOI.plotVolumeLayerFaces(geoV), view(3)
VOI.plotVolumeBoundaries(geoV), view(2)
%%
close all
%% Build the layered unstructured VOI grid
% The unstructured VOI grid includes a 2D well region (WR). The WR is 
% composed of a Cartesian region and two half-radial regions in xy plane, 
% which is used to connect the HW grid. 
% Generate the grid nodes for the 2D WR region. For the Cartesian region, 
% the X axis extends along the well trajectory:
%     ---------> X
%    |    -----------------------------------
%  Y |    --------- Well trajectory ---------
%    V    -----------------------------------
%
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

% Next, build the layered unstructured VOI grid. The refinement is allowed
% for each CPG layers.
% Define the number of refined layers for each VOI layer. Each of the four 
% VOI layers will be refined into two layers.
layerRf = [2, 2, 2, 2];

% The open-source triangle generator 'DistMesh' (Per-Olof Persson) is used
% to obtain high-quality triangles. 
% The scaled edge length function is defined as:
% h(p) = max(multiplier*d(p) +lIB, lOB)
% to let the point density increases from inner boundary to outer boundary
% lIB: average length of the inner boundary (outer-boundary of WR subgrid)
% lOB: average length of the outer boundary (clipped VOI boundary)

% Other parameters:
% 'maxIter' : The maximum number of distmesh iterations
% 'gridType': 'Voronoi' (default) | 'triangular'

% Reconstruct the CPG in VOI to layered unstructured grid
GV = VOI.ReConstructToUnstructuredGrid(WR, layerRf, ...
    'multiplier', 0.2, 'maxIter', 500, 'gridType', 'Voronoi');

% Show the VOI grid. We can plot the specified layers / surfaces by calling 
% 'GV.cells.layers==layer'/ 'GV.faces.surfaces==surface'
figure, axis off
arrayfun(@(layer)plotGrid(GV, GV.cells.layers==layer, ...
    'facecolor', rand(3,1)), 1:GV.layers.num)
view([-41, 62])

figure, axis off
arrayfun(@(surface)plotFaces(GV, GV.faces.surfaces==surface, ...
    'facecolor', rand(3,1)), 1:GV.layers.num+1)
view([-41, 62])
%%
close all
%% Build the layered radial HW grid
% The HW grid is built inside the Cartesian region of VOI grid. The logical
% indices of HW region should be specified.
%      1   ymin                     ymax   ny 
%    ----- ----- ----- ----- ----- ----- -----
%   |     |     |     |     |     |     |     |    1
%    ----- ----- ----- ----- ----- ----- -----
%   |     |  *  |  *  |  *  |  *  |  *  |     |   zmin
%    ----- ----- ----- ----- ----- ----- -----
%   |     |  *  |  *  |  *  |  *  |  *  |     |   
%    ----- ----- ----- ----- ----- ----- -----
%   |     |  *  |  *  |  *  |  *  |  *  |     |   zmax
%    ----- ----- ----- ----- ----- ----- -----
%   |     |     |     |     |     |     |     |    nz
%    ----- ----- ----- ----- ----- ----- -----
%   * = HW region
%   Remarks:    1 < ymin < ymax < ny (GV.surfGrid.cartDims(2))
%               1 < zmin < zmax < nz (GV.layers.num)
%          
%               [ymin, ymax, zmin, zmax]
regionIndices = [   2,    5,    2,    7];

% Define the HW region according to the GV, well, and regionIndices
HW = HorWellRegion(GV, well, regionIndices);

% Visualize the HW region
HW.showWellRegionInVOIGrid('showWellRgionCells', true);
view([-53, 15]), axis off

% Get the geometrical information of GV in HW region, including cells, 
% faces, boundary faces, nodes, boundary nodes, etc.
geoW = HW.allInfoOfRegion();
HW.plotRegionCells(geoW), axis equal, view([-85, 9])
HW.plotRegionLayerFaces(geoW), axis equal, view([-85, 9])

% Next, build the layered radial HW grid. We provide two types of grid 
% lines: 
% Type 1: 'pureCircular'  - The radial grid lines are pure circular
% % Parameter 'maxRadius' - Max radius of the radial grid
% % Parameter 'nRadCells' - Number of radial cells
radPara1 = struct(...
    'gridType'  , 'pureCircular', ...
    'maxRadius' , 2, ...
    'nRadCells' , 8);

% Type 2: 'gradual' - The radial grid lines vary from the circular line to 
%                     the rectangular line of a specified box gradually
% % Parameter 'boxRatio'  - Size ratio of the rectangular box to the outer 
% %                         boundary, [yRatio, zRatio]
% % Parameter 'nRadCells' - Number of radial cells, 2x1 double, 
% %                         [inbox, outbox]            
% % Parameter 'pDMult'    - Multiplier of pD of the outer-most angular line
%                           to pD of wellbore line. The outer-most line
%                           will be closer to the box with increasing pDMult
% % Parameter 'offCenter' - Whether the well is off-center in the radial 
% %                         grid
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

figure, axis equal tight off, view([-85, 9])
arrayfun(@(surface)plotFaces(GW, GW.faces.surfaces==surface, ...
    'facecolor', rand(3,1)), 1:GW.layers.num+1)
%%
close all
%% Collect the simulation data
% The standard ad-blackoil simulator requires 'G', 'rock', 'fluid', 
% 'model', 'schedule', and 'initState'. They should be redefined (except 
% the 'fluid') to conform with the global hybrid grid. The class 
% 'NearWellboreModel' is provided to collect the simulation data from 
% deck-format input.

% Define the NearWellboreModel by three subgrids, input deck and well
NWM = NearWellboreModel({GC, GV, GW}, deck, well);
%% Get the global hybrid grid
G = NWM.validateGlobalGrid();

% Show the global grid. We can plot the specified subgrids by calling 
% 'G.cells.grdID==i'
figure, hold on, axis off
plotGrid(G, G.cells.grdID == 1, 'facecolor', 'none')
plotGrid(G, G.cells.grdID == 2, 'facecolor', 'y')
view([-28, 20])

figure, hold on, axis off
plotGrid(G, G.cells.grdID == 2, 'facecolor', 'none')
plotGrid(G, G.cells.grdID == 3, 'facecolor', 'g')
view([-76, 70])

% Also, use 'G.cells.grdID == i & G.cells.layers == j' to plot layer j of
% subgrid i
figure, hold on, axis off
plotCellData(G, G.cells.centroids(:,1), ...
    G.cells.layers == 4 & G.cells.grdID == 1)
plotCellData(G, G.cells.centroids(:,1), ...
    G.cells.layers == 1 & G.cells.grdID == 2)
view([-3, 56])
%%
close all
%% Setup the fuild
fluid = NWM.setupFluid();
%% Get the rocks
% The rocks of three subgrids:
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

rockC = NWM.getCPGRockFromDeck();
rockV = NWM.getVOIRocksByInterp();

% View the permeability
figure
subplot(1,2,1), axis equal off
plotCellData(GC, rockC.perm(:,1), geoV.cells{1})
subplot(1,2,2), axis equal off
plotCellData(GV, rockV.perm(:,1), GV.cells.layers==1)

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
rock = NWM.getGlobalRock({rockC, rockV, rockW});

% View the permeability of global grid
figure, hold on, axis off
plotCellData(G, rock.perm(:,1), G.cells.layers == 4 & G.cells.grdID == 1)
plotCellData(G, rock.perm(:,1), G.cells.layers == 1 & G.cells.grdID == 2)
view([-3, 56])

figure, axis equal tight off, view([-85, 9])
plotCellData(G, rock.perm(:,1),G.cells.grdID == 3)
%%
close all
%% Get the transmissbility and neighbors
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
% * The T_nnc (transmissbility of non-neighbor connection) is used to
%   connect the three subgrids.
% * Updated GC and GV are obtained by 'removeCells'.
%
% Compute the transmissbility of the global grid T:
% T = [TC; TV; TW], corresponding to G.faces.neighbors
T = NWM.getTransGloGrid(rock);

% Generate the NNC
% First, compute the intersection relations between subgrids
intXn = NWM.computeBoundaryIntxnRelation();

% View intersection relation of the non-matching face
NWM.plotNonMatchingIntxnRelation(intXn ,15171)
% View intersection relation of the matching face
NWM.plotMatchingIntxnRelation(intXn ,2689)

% Generate the NNC
nnc = NWM.generateNonNeighborConn(intXn, rock, T);

% Get the assembled transmissbility and neighbors
% T_all = [T;                 nnc.T];
% N_all = [G.faces.neighbors; nnc.cells];
[T_all, N_all] = NWM.assembleTransNeighbors(T, nnc);
%%
close all
%% Setup simultaion model
% We use the 'GenericBlackOilModel' as the simulation model
model = NWM.setupSimModel(rock, T_all, N_all);
%% Get the simulation schedule
% Get the schedule from the production/injection control data from deck.
% Some fields in the well struct of the HW will be redefined, e.g. cells,
% 'WI'
schedule = NWM.getSimSchedule(model, 'refDepthFrom', 'deck');

% Show the well
figure, axis equal tight off, view([-85, 9])
plotGrid(G, G.cells.grdID == 3, 'facecolor', 'none')
plotWell(G, schedule.control(1).W(1))
%% Get the initial state by equilibrium initialization
initState = NWM.getInitState(model);
%% The packed output
% % Call 'packedSimData'to obtain all necessary simulation data of 
% % the near-wellbore model. The input is the rock of HW grid.
% [G, rock, fluid, model, schedule, initState] = NWM.packedSimData(rockW, 'refDepthFrom', 'deck');
% % Also, call 'packedCPGSimData' to obtain the simulation data of CPG.
% [GC, rockC, fluid, modelC, scheduleC, initStateC] = NWM.packedCPGSimData();
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