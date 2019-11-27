%% Near-wellbore modeling (NWM) workflow example
%
% This example aims to show the complete workflow for the NWM method. This
% includes building the volume of interest (VOI) grid, building the 
% horizontal well (HW) grid, and generating variables passed to the AD 
% simulators. Details about the grids and variable generation can be found
% in example 'nearWellBoreModelingGrids' and 'nearWellBoreModelingSim',
% respectively.

mrstModule add nwm ad-core ad-blackoil ad-props mrst-gui deckformat wellpaths upr

%% Read the ECLIPSE input deck
fn = fullfile(pwd, '..', 'data', 'NWM.data');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

%% Build the background Corner-point grid (CPG)
GC = initEclipseGrid(deck);
GC = computeGeometry(GC);

%% Define basic information of horizontal well (HW)
% Load well trajectory
fn = fullfile(pwd, '..', 'data', 'trajectory.mat');
load(fn)
% Number of well segments
ns = size(pW,1)-1;
% Define the well struct. 
well = struct(...
    'name'         , 'PROD', ...
    'trajectory'   , pW, ...
    'segmentNum'   , ns, ...
    'radius'       , 0.15 * ones(ns+1,1), ...
    'skinFactor'   , zeros(ns,1), ...
    'openedSegs'   , (1:ns));

%% Build the layered unstructured VOI grid
% Define the 2D boundary of VOI
pbdy = [240,   50;...
        160,   80;...
        120,  160;...
        150,  205;...
        230,  170;...
        280,   90];
    
% Define extra layers
nextra = [1, 1];

% Define the VOI according to the CPG, well, boundary and extra layers
VOI = VolumeOfInterest(GC, well, pbdy, nextra);

% Define parameters for the 2D well region grid 
WR.ny = 6; 
WR.ly = 12;
WR.na = 5;

% Define the number of refined layers for each VOI layer
layerRf = [2, 2, 2, 2];

% Reconstruct the CPG in VOI to layered unstructured grid
GV = VOI.ReConstructToUnstructuredGrid(WR, layerRf, ...
    'multiplier', 0.2, 'maxIter', 500, 'gridType', 'Voronoi');

%% Build the layered radial HW grid
% Define the logical indices of HW region
regionIndices = [   2,    5,    2,    7];

% Define the HW region according to the GV, well, and regionIndices
HW = HorWellRegion(GV, well, regionIndices);

% Define the parameters for building the radial grid
radPara2 = struct(...
    'gridType'  , 'gradual', ...
    'boxRatio'  , [0.6, 0.6], ...
    'nRadCells' , [7, 2], ...
    'pDMult'    , 10, ...
    'offCenter' , true);

% Reconstruct the VOI grid in HW region to layered radial grid
GW = HW.ReConstructToRadialGrid(radPara2);

%% Generate variables for the AD simulator
% Define the 'NearWellboreModel' by three subgrids, input deck and well
NWM = NearWellboreModel({GC, GV, GW}, deck, well);

% Define a simple rock for HW grid
nclayer = GW.cells.num / GW.layers.num;
permW = linspace(400, 500, GW.layers.num) * (milli*darcy);
permW = repmat(permW, nclayer, 1);
poroW = linspace(0.18, 0.2, GW.layers.num);
poroW = repmat(poroW, nclayer, 1);
rockW.perm = [permW(:), permW(:), permW(:)];
rockW.poro = poroW(:);

% Get the variables for the NWM grid
[G, rock, fluid, model, schedule, initState] = NWM.packedSimData(rockW);

% Get the variables for the CPG
[GC, rockC, ~, modelC, scheduleC, initStateC] = NWM.packedCPGSimData();

%% Run the simulator
[wellSols, states, report] = simulateScheduleAD(initState, model, schedule);

