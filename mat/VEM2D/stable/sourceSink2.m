clc; clear all; close all;

addpath('../');
addpath('../../');
addpath('~/Documents/master/pebiGridding/voronoi2D/')
run('../../startup.m');

n = 20;

G = cartGrid([n,n]);
G = computeGeometry(G);
G = mrstGridWithFullMappings(G);
G = sortEdges(G);


faultFaces = 1:G.faces.num;
faultFaces = faultFaces(G.faces.centroids(:,1) == ceil(n/2) ...
                        & G.faces.centroids(:,2) < n*3/4  ...
                        & G.faces.centroids(:,2) > n/4);
G = VEM2D_makeInternalBoundary(G, faultFaces);
f = zeros(G.cells.num,1);

G = computeVEM2DGeometry(G,f,1,1);

sourceCoords = [n*.25,n/2];
source = sum(bsxfun(@minus, G.cells.centroids, sourceCoords).^2,2);
source = find(source == min(source));

sinkCoords = [.6,.6];
sink = sum(bsxfun(@minus, G.cells.centroids, sinkCoords).^2,2);
sink = find(sink == min(sink));

f = @(X) zeros(size(X,1),1);

Q = 10;
src = addSource([], source, Q);
% src = addSource(src, sink, -Q);

boundaryEdges = find(G.faces.neighbors(:,1) == 0 | G.faces.neighbors(:,2) == 0);
isExternal = G.faces.centroids(boundaryEdges,1) == 0 | ...
        G.faces.centroids(boundaryEdges,1) == n | ...
        G.faces.centroids(boundaryEdges,2) == 0 | ...
        G.faces.centroids(boundaryEdges,2) == n;
isInternal = ~isExternal;

bc = VEM_addBC(G, [], boundaryEdges(isExternal), 'pressure', 0);
bc = VEM_addBC(G, bc, boundaryEdges(isInternal), 'flux', 0);
            
sol1 = VEM2D_v3(G,f,1,bc,'src', src);
sol1 = cellAverages(G,sol1);

% plotVEM(G, sol1.nodeValues, 'dof')

sol2 = VEM2D_v3(G,f,2,bc, 'src', src);

plotVEM(G,[sol2.nodeValues; sol2.edgeValues; sol2.cellMoments], '')

figure;
plotCellData(G, sol1.cellAverages);
plotFaces(G,boundaryEdges(isInternal),'EdgeColor', 'r')
colorbar;
figure;
plotCellData(G, sol2.cellMoments)
plotFaces(G,boundaryEdges(isInternal),'EdgeColor', 'r')
colorbar;
