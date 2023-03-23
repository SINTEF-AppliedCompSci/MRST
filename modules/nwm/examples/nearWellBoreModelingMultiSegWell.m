%% Coupled model of the near-wellbore model (NWM) and multi-segment well (MSW)
%
% This example displays the coupled model of NWM and MSW, which aims to 
% describe both the near-wellbore flow and flow inside the wellbore. The
% MSW requires the node and segment definitions, which is similar to the
% cell and face definitions if we consider the MSW running on a virtual 
% grid domain. Therefore, the 'MultiSegWellNWM' will automatically generate
% the node and segment definitions through building a 'wellbore grid'. The 
% geometrical information was obtained by grid geometries, e.g. 
% 'cells.centroids', 'cells.volumes', and the topology of the segments 
% corresponds to 'faces.neighbors'. The pressure drop calculation adopts 
% the mrst built-in wellbore friction model.

clear
mrstModule add nwm ad-core ad-blackoil ad-props mrst-gui deckformat wellpaths upr

%% Read the ECLIPSE input deck
fn = fullfile('data', 'MSW.data');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

%% Build the background Cartesian grid
% Here, we build a Cartesian grid as the background grid
GC = initEclipseGrid(deck);
GC = computeGeometry(GC);

%% Define basic information of horizontal well (HW)
% Define the well trajectory from deck
compi = cell2mat( deck.SCHEDULE.control.COMPDAT(:,2) );
compj = cell2mat( deck.SCHEDULE.control.COMPDAT(:,3) );
compk = cell2mat( deck.SCHEDULE.control.COMPDAT(:,4) );
[dx, dy, dz]   = deal(unique(deck.GRID.DX), unique(deck.GRID.DY), unique(deck.GRID.DZ));
tops = deck.GRID.TOPS(1);

pWx = dx * linspace(min(compi)-1, max(compi), 15)';
pWy = dy * ( unique(compj) - 0.5 ) * ones( size(pWx) );
pWz = dz * ( unique(compk) - 0.5 ) * ones( size(pWx) ) + tops;
pW = [pWx, pWy, pWz];
% Number of well segments
ns = size(pW,1)-1;
% Define the well structure
well = struct(...
    'name'         , 'PROD', ...
    'trajectory'   , pW, ...
    'segmentNum'   , ns, ...
    'radius'       , 0.1 * ones(ns+1,1), ...
    'skinFactor'   , zeros(ns,1), ...
    'openedSegs'   , 1:ns, ...
    'isMS'         , true, ...
    'roughness'    , 1e-2 * ones(ns,1));
% We set the roughness and the well production rate to quite large values 
% to show the difference between MSW and simple well.

%% Build the layered unstructured VOI grid
% Define the 2D boundary of VOI
pbdy = 100* [ 4,   10;...
             21,   10;...
             21,   15;...
              4,   15];
          
% Define extra layers
nextra = [1, 1];

% Define the VOI according to the CPG, well, boundary and extra layers
VOI = VolumeOfInterest(GC, well, pbdy, nextra);

% Define parameters for the 2D well region grid 
VOI.maxWellSegLength2D()
WR = struct('ly', [15, 80], 'ny', [6, 4], 'na', 5);

% Define the number of refinement layers for each VOI layer
VOI.volumeLayerNumber()
layerRf = [2, 2, 2];

% Reconstruct the CPG in VOI to layered unstructured grid
GV = VOI.ReConstructToUnstructuredGrid(WR, layerRf, ...
    'multiplier', 0.2, 'maxIter', 500, 'gridType', 'Voronoi');

% Plot the VOI grid
figure, axis equal tight off, view(3)
arrayfun(@(layer)plotGrid(GV, GV.cells.layers==layer, ...
    'facecolor', rand(3,1)), 1:GV.layers.num)
title('Layers of the VOI grid')

%% Build the layered radial HW grid
% Define the logical indices of HW region
regionIndices = [   3,    8,    2,    5];

% Define the HW region according to GV, well, and regionIndices
HW = HorWellRegion(GV, well, regionIndices);

HW.showWellRegionInVOIGrid('showWellRgionCells', true);
view([-53, 15]), axis equal off
title('The HW region in VOI grid')


% Define the parameters for building the radial grid
radPara = struct(...
    'gridType'  , 'gradual', ...
    'boxRatio'  , [0.6, 0.6], ...
    'nRadCells' , [10, 2], ...
    'pDMult'    , 15, ...
    'offCenter' , true);

% Reconstruct the VOI grid in HW region to layered radial grid
GW = HW.ReConstructToRadialGrid(radPara);

% Plot the HW grid
figure, axis equal tight off, view(3)
arrayfun(@(layer)plotGrid(GW, GW.cells.layers==layer, ...
    'facecolor', rand(3,1)), 1:GW.layers.num)
title('Layers of the HW grid')

%% Collect the simulation data
% Define the Multi-Segment well NWM by three subgrids, input deck and well
MSW = MultiSegWellNWM({GC, GV, GW}, deck, well);

% Get the global hybrid grid
G = MSW.validateGlobalGrid();

% Initialize the AD fluid
fluid = MSW.setupFluid();

% Make the subrocks
rockC = MSW.getCPGRockFromDeck();  % CPG rock
rockV = MSW.getVOIRocksByInterp(); % VOI grid rock
% Define a homogeneous rock for HW grid
rockW.perm = ones(GW.cells.num, 3) * rockC.perm(1,1);
rockW.poro = ones(GW.cells.num, 1) * 0.2;

% Get the rock for global grid
rock = MSW.getGlobalRock({rockC, rockV, rockW});

% Compute the transmissbility and neighborship
T = MSW.getTransGloGrid(rock);
intXn = MSW.computeIntxnRelation();
nnc = MSW.generateNonNeighborConn(intXn, rock, T);
[T_all, N_all] = MSW.assembleTransNeighbors(T, nnc);

% Setup simulation model
% Note the MSW model now only supports the 'ThreePhaseBlackOilModel'
model = MSW.setupSimModel(rock, T_all, N_all);

%% Get the initial state by equilibrium initialization
initState = MSW.getInitState(model);

%%
% Above procedures (gridding of VOI and HW, collections of G, fluid, rock, 
% model, and initState) are same with the basic near-wellbore model

%% Convert the simulation schedule
% We build a 1D 'wellbore grid' in the void wellbore space and it conforms
% with the reservoir grid. The wellbore grid has 'ns' cells (nodes) and 
% 'ns-1' internal faces (segments). 
% View the wellbore grid
gW  = MSW.wellboreGrid;
wc  = MSW.getWellCells();
figure, hold on, axis tight off
plotGrid(gW)
plotGrid(G, wc, 'facecolor', 'none')
legend('Wellbore grid', 'Well cells of reservoir grid')
view(3)

% The cells/internal faces are equivalent to the nodes/segments in the
% multi-segment well model. Each node accepts inflow from reservoir cells
% in corresponding grid layer. The segment connects its forward node and
% backward node:
%   
%     o    o     o     o     o     o     o     o     o   reservoir cells
%     |    |     |     |     |     |     |     |     |   inflow to nodes
%     * -- *  -- *  -- *  -- *  -- *  -- *  -- *  -- *   nodes and segments
%     |    |     |     |     |     |     |     |     |
%     o    o     o     o     o     o     o     o     o
%
% Specially, the nodes are generated as:
%  'depth'       : gW.cells.centroids(:,3)
%  'vol'         : gW.cells.volumes
%  'cell2node'   : From grid layer indices
% The segments are generated as:
%  'topo'        : gW.faces.neighbors(internal, :)
%  'diam'        : From well info structure
%  'roughness'   : From well info structure
%  'length'      : Length between forward node and backward node 
%  'flowModel'   : mrst 'wellBoreFriction' model
nodes    = MSW.generateNodes()
segments = MSW.generateSegments()
% Plot the nodes and reservoir cells associated with segment 1:2
MSW.plotSegments(nodes, segments, 1:2)
axis tight off, view([-16, 45])
title('The segments, nodes and reservoir cells connected to nodes')

% Convert the schedule from CPG to global grid
% Note the reference depth has been set to the depth of the top node
% (i.e. first cell in wellbore grid)
scheduleMS = MSW.getSimSchedule(model);

% We can also get the schedule without multi-segment well definition
schedule = MSW.getSimSchedule(model, 'returnMS', false);

%% Run the simulations with and without multi-segment well
[wellSolsMS, statesMS, reportMS] = simulateScheduleAD(initState, model, scheduleMS);
[wellSols, states, report] = simulateScheduleAD(initState, model, schedule);

%% Compare the well solutions
plotWellSols({wellSols, wellSolsMS}, report.ReservoirTime)

%% Compare the liquid production along the well
nA = GW.radDims(1);
x = ( pW(1:end-1,1) + pW(2:end,1) ) / 2;
[xx, tt] = meshgrid(x,report.ReservoirTime/day);

[fluxMS, flux] = deal( zeros(numel(wellSols), ns) );
for i = 1 : numel(wellSols)
    fMS = -sum(wellSolsMS{i}.flux,2);
    fMS = reshape(fMS, nA, []);
    fluxMS(i, :) = sum(fMS, 1);
    f = -sum(wellSols{i}.flux,2);
    f = reshape(f, nA, []);
    flux(i, :) = sum(f, 1);
end

figure,hold on
surfWithOutline(xx, tt, fluxMS/(meter^3/day));
surfWithOutline(xx, tt, flux/(meter^3/day));
hold off, box on, axis tight
shading interp
set(gca,'Projection','Perspective');
hc = get(gca, 'Children');
set(hc(1:4), 'Color', 'r')
xlabel('Distance from HW heel (m)')
ylabel('Time (d)')
zlabel('Liquid production (m^3/day)')
view(3)

%% Plot the node pressures
ts = 20;
ws = wellSolsMS{ts};
nPres = [ws.bhp; ws.nodePressure];
L = nodes.coords(:,1);
L = L - L(1);
figure
plot(L, nPres/barsa, 's-')
xlabel('Distance from top node (m)')
ylabel('Node pressure (barsa)')

state = statesMS{end};
figure, axis tight off, view(3)
plotCellData(G, state.pressure/barsa, wc)
plotCellData(gW, nPres/barsa)
title('Pressures of the near-wellbore grid and wellbore grid')
