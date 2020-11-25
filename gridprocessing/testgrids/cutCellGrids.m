% Some simple tests
mrstModule add libgeometry

%% Simple stress-test interatively cutting a 1x1x1 box
%  input type 1: single coordinate + cut-plane normal
figure, axis off
G = addBoundingBoxFields(computeGeometry(cartGrid([1 1 1], [20, 20, 20])));
plotGrid(G,'FaceAlpha', .1, 'EdgeAlpha', .5), view(3), camlight
col = {'r', 'b', 'g', 'm', 'c', 'y'};
nIts = 50;
cc   = nan(nIts,1);
for k = 1:nIts
    [G, ix, g2] = sliceGrid(G, 10 + 18*rand(1,3)-9, 'normal', 20*rand(1,3)-10,...
        'radius', inf*(1 + 10*rand(1,3)), 'topoSplit', true);
    if ~isempty(ix)
        plotFaces(G, ix.new.faces==3, 'FaceColor', col{randi(6)}, 'FaceAlpha', 1);
    end
    title(sprintf('Cut no: %d, cells: %d', k, G.cells.num));
    chk = checkGrid(G);
    if ~chk
        return
    end
    drawnow
end
figure, plotCellData(G, G.cells.volumes), view(3), camlight

%% Cut along (piecewise linear) curve  
% input type 2: curve-coordinates + cutting direction
% limitation is that slice cannot enter and exit cell trough the same face, 
% and hence the curve needs to be somewhat smooth wrt to grid resolution 
% (this is not checked for, but will produce problems.)
figure, axis off 
G = addBoundingBoxFields(computeGeometry(cartGrid([20 20 5], [20, 20, 5])));
plotGrid(G,'FaceAlpha', .1, 'EdgeAlpha', .5), view(3), camlight
x = (-1:.5:21)';
points = [x, 10+5*cos(.2*x), .2*x];
cutDir = [0 -.5 1];
[G, ix, g21] = sliceGrid(G, points, 'cutDir', cutDir);
plotFaces(G, ix.new.faces==3, 'FaceColor', 'm', 'FaceAlpha', 1, 'EdgeAlpha', .2);
%
points = [x, x, 2.8 + 2*sin(.7*x) + .05*x];
cutDir = [1, -1, 0.1];
[G, ix, g22] = sliceGrid(G, points, 'cutDir', cutDir, 'radius', inf);
plotFaces(G, ix.new.faces==3, 'FaceColor', 'c', 'FaceAlpha', 1, 'EdgeAlpha', .2);
chk = checkGrid(G)
%%
% 2D grid-slices are not yet connected to single grid
mrstModule add incomp
g22 = repairNormals(computeGeometry(g22));
rock = makeRock(g22, 100*milli*darcy, 1);
fluid = initSingleFluid('mu',1*centi*poise, 'rho', 1014*kilogram/meter^3);
[~, f1] = min(g22.faces.centroids(:,1:2)*[1 1]');
[~, f2] = max(g22.faces.centroids(:,1:2)*[1 1]');
bc = addBC([], [f1;f2], 'pressure', [500; 50]*barsa);
rs = incompTPFA(initState(g22,[],0), g22, ones(nnz(g22.faces.neighbors>0),1)*darcy, fluid, 'bc', bc);
figure, axis off 
plotCellData(g22, rs.pressure), view([1 -2 10]), camlight


%% Slicing of 2D grids follow same procedure (needs 3D coordinates)
G = cartGrid([50 50], [100, 100]);
G.nodes.coords(:,3) = 0;
G = repairNormals(computeGeometry(G));
G = repairNormals(G);
figure, 
for k = 1:50
    [G, ix] = sliceGrid(G, [100*(rand(1,2)), 0],  'normal', [rand(1,2)-.5, 0],...
        'radius', (1 + 50*rand));
    if ~isempty(ix)
        cla, plotCellData(G, G.cells.volumes, 'EdgeAlpha', .2);
        drawnow
    end
end
    
        
