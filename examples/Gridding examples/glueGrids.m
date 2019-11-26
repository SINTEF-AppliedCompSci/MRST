%% Demonstration of gluing the Voronoi grid to Cartesian grid
%
% This example demonstrates how to glue the Voronoi grid to the Cartesian 
% grid along a complex boundary clipped by the Cartesian edges (faces).
% This example, also most of the 'near-wellbore modeling' 2D gridding 
% process, adopt the constructor 'tessellationGrid' from mrst grid factory. 
% The 'tessellationGrid' accepts the combination of points (geometrical 
% input) and connectivity list (topological inputs) and returns a complete
% grid structure. The connectivity list provides the indices of nodes
% specifing each polygonal cell, i.e. the 'tessellationGrid' restores the 
% complete grid-object sequence of nodes->faces->cells from nodes->cells.
% The connections of two grids is then quite easy: assemble the two point 
% sets and connectivity lists to obtain the global ones. However, the
% common points from the two point sets should be merged and the 
% connectivity list should be mapped, otherwise the flow interactions 
% between the two grids are unavailable. An good alternative to merging the
% points is generating the 'non-neighbor connection (NNC)' that is quite
% practical to connect the non-neighboring cells. This example displays
% both the two grid gluing methods and compares the simulation results.
% For another, the example shows how to generate the Voronoi grid under the
% fix-point constraint and how to clip the infinte Voronoi diagram.
%
% Load necessary modules
mrstModule add ad-core ad-blackoil ad-props
%% Get boundary information
% Build the Cartesian grid
GC = cartGrid([25, 25], [200, 200]);
GC = computeGeometry(GC);

% Define the well trajectory and volume of interest (VOI) boundary polygon
% The well should be located inside the VOI boundary.
ns = 12;
th = linspace(0.85, 0.7, ns+1)' * pi;
traj = [150*cos(th)+200, 150*sin(th)];

pbdy = [136, 150;
        145,  95;
         90,  30;
         50,  50;
         45, 105;
         90, 160];

% Get cells inside and outsied the VOI
cCtro = GC.cells.centroids;
in = inpolygon(cCtro(:,1), cCtro(:,2), pbdy(:,1), pbdy(:,2));
cI = find( in );
cO = find(~in );

% Get the sorted boundary nodes and boundary cells (in 'counterclockwise').
% The actual VOI boundary is aligned to the grid edges (faces) and 
% therefore lead the nonconvex polygons. So, we first sort the convex 
% polygon specified by the centroids of boundary edges, and then obtain the 
% sorted boundary nodes and cells from associated faces.
[bnv, bcv] = demoGetBdyNodesCells(GC, cI);

% Plot the VOI boundary
figure, hold on, axis equal tight off
plotGrid(GC, cO, 'facecolor', 'none')
plotGrid(GC, cI, 'facecolor', 'y')
demoPlotLine(traj, 'ko-', 'b', 4)
demoPlotPoly(pbdy, 'k^-', 'r', 5)
demoPlotPoly(GC.nodes.coords(bnv,:), 'ks-', 'g', 4)
legend('GC outside the VOI', 'GC inside the VOI', 'Well path' ,...
    'Specified VOI boundary','Clipped VOI boundary')
%% Get the well-region points and connectivity list
% The WR is composed of a Cartesian region and two half-radial regions.
ly = 12;
ny = 8;
na = 6;
pw0 = arrayfun(@(ii)pointsSingleWellNode(traj, ly, ny, na, ii), (1:ns+1)');
pw  = [vertcat(pw0.cart); vertcat(pw0.rad)];
% Make sure all well-region points are loacted inside the VOI
assert( all(inpolygon(pw(:,1), pw(:,2), GC.nodes.coords(bnv,1), ...
    GC.nodes.coords(bnv,2))) );
[tw, ~, bnw] = getConnListAndBdyNodeWR2D(pw0, ny, na);

figure, hold on, axis equal tight off
plotGrid(GC, cO, 'facecolor', 'none')
demoPlotPoly(GC.nodes.coords(bnv, :), 'ks-', 'g', 4)
demoPlotPoly(pw(bnw, :), 'ko-', 'y', 3)
legend('GC outside the VOI', 'GWR (WR grid)', 'VOI boundary' ,'WR boundary')
%%
pob  = GC.nodes.coords(bnv,:);
pob2 = GC.cells.centroids(bcv, :);
pib  = pw(bnw, :);
[pIn, pOut, R] = demoComputeAuxPts(pw, bnw, 0.35);
pib2 = pOut;

figure, hold on, axis equal tight off
demoPlotPoly(pob,  'ks-', 'g', 4)
demoPlotPoly(pib,  'ko-', 'y', 3)
demoPlotPoly(pob2, 'ks-', 'r', 4)
demoPlotPoly(pib2, 'ko-', 'b', 3)
legend('VOI boundary' ,'WR boundary', 'VOI boundary for tri. pts. generation' ...
    ,'WR boundary for tri. pts. generation')
%%
% Generate basic points from distmesh
[pdis, tdis] = passToDistmesh(pib2, pob2, 0.2, 500, 'pIBRadius', R);

pauxw = pIn;
pauxv = GC.cells.centroids(cO, :);

% Get Voronoi points and connectivity list
pall = [pdis; pauxv; pauxw];
[pVor, tVor] = voronoin(pall, {'Qbb','Qz'});
%%
figure, hold on, axis equal off
voronoi(pall(:,1), pall(:,2), '-')
demoPlotPoly(pdis,  'ro',  'r', 2)
demoPlotPoly(pauxv, 'bo',  'b', 2)
demoPlotPoly(pauxw, 'k^',  'k', 3)
demoPlotPoly(pob,   'ks-', 'g', 4)
demoPlotPoly(pib,   'ko-', 'y', 3)
legend('', 'Initial Voronoi digram', 'Basic sites' ,'Auxiliary sites (CO)', ...
    'Auxiliary sites (WR)',  'VOI boundary' ,'WR boundary')
xlim([-20, 220])
ylim([-20, 220])
%%
[pv, tv, bnv2] = demoHandleVoronoiDiag(pVor, tVor, pib, pob, pw, tw, bnw);

tv = sortPtsCounterClockWise(pv, tv);
GV = tessellationGrid(pv, tv);
GV = computeGeometry(GV);

figure, hold on, axis equal tight off
plotGrid(GC, cO, 'facecolor', [.0, .44, .74])
plotGrid(GV, 'facecolor', [.85, .32, .09])
legend('GC outside the VOI', 'GV (New VOI grid)')

