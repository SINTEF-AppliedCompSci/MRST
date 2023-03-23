%% The grids in the Near-wellbore modeling (NWM) method
% 
% Typically, the horizontal well (HW) is built in the Corner-point grid 
% (CPG) or Cartesian grid by the standard well model (Peaceman's model) 
% that inserts the virtual well into coarse reservoir grid. However, the 
% standard model gives a grid-based well trajectory, which will be rather 
% unrigorous when the real well trajectory fails to follow the grid 
% orientation. For another, some field operations require high-resolution 
% flow descriptions in the well vicinity. The coarse well cells always 
% unable to provide such descriptions due to the linear flow approximation 
% and low-resolution rock properties. To this end, this example 
% demonstrates the local grid reconstruction in the near-wellbore region of
% horizontal well. 
%
% The global grid will consist of three subgrids that adopt different
% gridding strategies:
% -------------------------------------------------------------------------
% | Grid | Description | Type                 | Constructor               |
% |-----------------------------------------------------------------------|
% | GC   | Background  | Corner-point or      | initEclipseGrid           |
% |      | grid        | Cartesian            |                           |
% |-----------------------------------------------------------------------|
% | GV   | VOI grid    | Unstructured,        | tessellationGrid +        |
% |      |             | Vertically layered   | makeLayeredGridNWM        |
% |-----------------------------------------------------------------------|
% | GW   | HW grid     | Structured, Radial,  | buildRadialGrid +         |
% |      |             | Horizontally layered | makeLayeredGridNWM        |
% -------------------------------------------------------------------------
%
% Remarks:
% * This module now only supports the modeling of horizontal wells. The
%   modeling of vertical wells is under development.
% * The volume of interest (VOI) cannot cover the faults.
% * The package 'distmesh' comes from module 'upr'. Note module 'hfm' also
%   includes the 'distmesh', but the iterations will be slow due to the  
%   pure MATLAB version of the DSEGMENT routine introduced in 'hfm'.

clear
mrstModule add nwm deckformat wellpaths upr

%% Read the ECLIPSE input deck
fn = fullfile('data', 'NWM.data');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

%% Build the background Corner-point grid (CPG)
GC = initEclipseGrid(deck);
GC = computeGeometry(GC);

%% Define basic information of horizontal well (HW)
% The well trajectory is specified by a set of discrete 3D well points 
% (in xyz format) which divides the HW into multiple segments.
% Load well trajectory
fn = fullfile('data', 'trajectory.mat');
load(fn)
% Number of well segments
ns = size(pW,1)-1;

% Define the well structure 
%  'name':         Well name, should match the well name list in deck 
%  'trajectory':   Well trajectory, 3D points in xyz format
%  'segmentNum':   Number of well segment, equals to n_wellpoints - 1
%  'radius':       Wellbore (Casing) radius per well segment
%  'skinFactor':   Skin factor per well point
%  'openedSegs':   Opened segments (allow fluid to flow into the wellbore)
% The coupling of the NWM with multi-segment well requires: 
%  'isMS'      :   Indicating the multi-segment well definition
%  'roughness' :   Roughness per well segment
% (See 'nearWellBoreModelingMultiSegWell')
well = struct(...
    'name'         , 'PROD', ...
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
% nextra(1): Number of layers above the HW layers
% nextra(2): Number of layers below the HW layers
nextra = [1, 1];

% Define the VOI according to the CPG, well, boundary and extra layers
VOI = VolumeOfInterest(GC, well, pbdy, nextra);

% Get the geometrical information of CPG in VOI, including cells, 
% layer-faces, boundary faces, nodes, boundary nodes, etc.
geoV = VOI.allInfoOfVolume();

% Show the VOI boundary, cells, and faces.
VOI.plotVolumeCells(geoV), view(3)
VOI.plotVolumeLayerFaces(geoV), view(3)
VOI.plotVolumeBoundaries(geoV), view(2)

%% Build the layered unstructured VOI grid
% The unstructured VOI grid includes a 2D well region (WR). The WR is 
% composed of a Cartesian region and two half-radial regions in xy plane, 
% which is used to connect the HW grid. For the Cartesian region, the X 
% axis extends along the well trajectory:
%     ---------> X
%    |    -----------------------------------
%  Y |    --------- Well trajectory ---------
%    V    -----------------------------------

% ly: The size of Cartesian region in Y direction, better to be larger than
% the well-segment length
% ny: The number of Cartesian cells in Y direction
% na: The number of angular cells in radial region
VOI.maxWellSegLength2D()

% if numel(ly) == 1 && numel(ny) == 1 :
% Uniform  distribution
WR = struct('ly', 15, 'ny', 10, 'na', 5);
VOI.plot2DWRSubGrid(WR)

% if numel(ly) == 2 && numel(ny) == 2 :
% ly(1) ny(1) : uniform  distribution in inner region
% ly(2) ny(2) : logarithmic distribution in outer region
WR = struct('ly', [5,10], 'ny', [6,4], 'na', 5);
VOI.plot2DWRSubGrid(WR)

% Next, build the layered unstructured VOI grid. The refinement is allowed
% for each CPG layer.
% Define the number of refinement layers for each VOI layer. The dimension 
% of 'layerRf' should be equal to the number of VOI layers.
VOI.volumeLayerNumber()
% Each of the four VOI layers will be refined into two sublayers.
layerRf = [2, 2, 2, 2];

% The open-source triangle generator 'DistMesh' (Per-Olof Persson) is used
% to obtain high-quality triangles. The scaled edge length function is 
% defined as:
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
figure, axis off, view([-41, 62])
arrayfun(@(layer)plotGrid(GV, GV.cells.layers==layer, ...
    'facecolor', rand(3,1)), 1:GV.layers.num)
title('Layers of the VOI grid')

figure, axis off, view([-41, 62])
arrayfun(@(surface)plotFaces(GV, GV.faces.surfaces==surface, ...
    'facecolor', rand(3,1)), 1:GV.layers.num+1)
title('Surfaces of the VOI grid')

%% Build the layered radial HW grid
% The HW grid is built inside the Cartesian region of VOI grid. The logical
% indices of HW region should be specified:
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
regionIndices = [   3,    8,    2,    7];

% Define the HW region according to GV, well, and regionIndices
HW = HorWellRegion(GV, well, regionIndices);

% Visualize the HW region
HW.showWellRegionInVOIGrid('showWellRgionCells', true);
view([-53, 15]), axis off
title('The HW region in VOI grid')

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
    'boxRatio'  , [0.6, 0.6], ...
    'nRadCells' , [7, 2], ...
    'pDMult'    , 10, ...
    'offCenter' , true);

% Reconstruct the VOI grid in HW region to layered radial grid
GW = HW.ReConstructToRadialGrid(radPara2);

% Visualize the raidal grid in HW region
HW.showWellRegionInVOIGrid('showWellRgionCells', false);
plotGrid(GW, GW.cells.layers==1, 'facecolor', rand(3,1))
view([-53, 15]), axis off
title('The reconstructed HW region in VOI grid')

% Show the HW grid
figure, axis equal tight off, view([-85, 9])
arrayfun(@(layer)plotGrid(GW, GW.cells.layers==layer, ...
    'facecolor', rand(3,1)), 1:GW.layers.num)
title('Layers of the HW grid')

figure, axis equal tight off, view([-85, 9])
arrayfun(@(surface)plotFaces(GW, GW.faces.surfaces==surface, ...
    'facecolor', rand(3,1)), 1:GW.layers.num+1)
title('Surfaces of the HW grid')