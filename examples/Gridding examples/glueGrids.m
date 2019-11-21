clear
close all
%%
GC = cartGrid([25, 25], [200, 200]);
GC = computeGeometry(GC);

ns = 12;
th = linspace(0.18, 0.32, ns+1)' * pi;
R = 150;
traj = [R*cos(th), R*sin(th)];

pbdy = [50,150; 60,80;120,36;160,60;150,110;100,170];

cCtro = GC.cells.centroids;
in = inpolygon(cCtro(:,1), cCtro(:,2), pbdy(:,1), pbdy(:,2));
CI = find( in );
CO = find(~in );

[bnv, bcv] = demoGetBdyNodesCells(GC, CI);

figure, hold on, axis equal tight off
plotGrid(GC, CO, 'facecolor', 'none')
plotGrid(GC, CI, 'facecolor', 'y')
demoPlotLine(traj,   'ko-', 'b', 4)
demoPlotPoly(pbdy, 'k^-', 'r', 5)
demoPlotPoly(GC.nodes.coords(bnv,:), 'ks-', 'g', 4)
legend('GC outside the VOI', 'GC inside the VOI', 'Well path' ,...
    'Specified VOI boundary','Clipped VOI boundary')
%%
ly = 12;
ny = 8;
na = 6;
pw0 = arrayfun(@(ii)pointsSingleWellNode(traj, ly, ny, na, ii), (1:ns+1)');
pw  = [vertcat(pw0.cart); vertcat(pw0.rad)];
assert(all(inpolygon(pw(:,1), pw(:,2), pbdy(:,1), pbdy(:,2))),...
    ['Points outside the boundary were detected, please reduce ',...
    'the size of Cartesian region']);
[tw, ~, bnw] = getConnListAndBdyNodeWR2D(pw0, ny, na);

GWR = tessellationGrid(pw, tw);

figure, hold on, axis equal tight off
plotGrid(GC, CO, 'facecolor', 'none')
plotGrid(GWR, 'facecolor', 'c')
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
demoPlotPoly(pob, 'ks-', 'g', 4)
demoPlotPoly(pib, 'ko-', 'y', 3)
demoPlotPoly(pob2, 'ks-', 'r', 4)
demoPlotPoly(pib2, 'ko-', 'b', 3)
legend('VOI boundary' ,'WR boundary', 'VOI boundary for tri. pts. generation' ...
    ,'WR boundary for tri. pts. generation')
%%
% Generate basic points from distmesh
[pdis, tdis] = passToDistmesh(pib2, pob2, 0.2, 500, 'pIBRadius', R);

pauxW = pIn;
pauxV = GC.cells.centroids(CO, :);

% Get Voronoi points and connectivity list
pall = [pdis; pauxV; pauxW];
[pVor, tVor] = voronoin(pall, {'Qbb','Qz'});
%%
figure, hold on, axis equal off
voronoi(pall(:,1), pall(:,2), '-')
demoPlotPoly(pdis,  'ro', 'r', 2)
demoPlotPoly(pauxV, 'bo', 'b', 2)
demoPlotPoly(pauxW, 'k^', 'k', 3)
demoPlotPoly(pob, 'ks-', 'g', 4)
demoPlotPoly(pib, 'ko-', 'y', 3)
legend('', 'Initial Voronoi digram', 'Basic sites' ,'Auxiliary sites (CO)', ...
    'Auxiliary sites (WR)',  'VOI boundary' ,'WR boundary')
xlim([-20, 220])
ylim([-20, 220])
%%
[p, t, bdyID] = demoHandleVoronoiDiag(pVor, tVor, pib, pob, pw, tw, bnw);

t = sortPtsCounterClockWise(p, t);
GV = tessellationGrid(p, t);

figure, hold on, axis equal tight off
plotGrid(GC, CO, 'facecolor', [.0, .44, .74])
plotGrid(GV, 'facecolor', [.85, .32, .09])
legend('GC outside the VOI', 'GV (New VOI grid)')
%%
[GC, ~, ~, mapn] = removeCells(GC, CI);
%%
bnv = arrayfun(@(n)find(mapn == n), bnv);
%%
nnv = GV.nodes.num;
%%
pc = GC.nodes.coords;
[cnodes, pos] = gridCellNodes(GC, (1:GC.cells.num));




cnodes = arrayfun(@(c)gridCellNodes(GC, c), (1:GC.cells.num)', 'UniformOutput' ,false);
cnodes = sortPtsCounterClockWise(pC, cnodes);
%%


