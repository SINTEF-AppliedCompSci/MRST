clc; clear all; close all;

% addpath('../');
% addpath('../../');
% addpath('~/Documents/master/pebiGridding/voronoi2D/')
% run('../../startup.m');

n = 3;

% G = cartGrid([n,n],[n,n]);



gridLim = [1,1];
G = cartGrid([n,n],gridLim);
% G = unitSquare(n,n);

G = computeGeometry(G);
G = sortEdges(G);
G = mrstGridWithFullMappings(G);


faultFaces = 1:G.faces.num;
faultFaces = faultFaces(min((abs(G.faces.centroids(:,1) - .5*gridLim(1)))) ...
                          == abs(G.faces.centroids(:,1)-.5*gridLim(1)) ...
                        & G.faces.centroids(:,2) < .9*gridLim(2)  ...
                        & G.faces.centroids(:,2) > .1*gridLim(2));

G = VEM2D_makeInternalBoundary(G, faultFaces);

f = zeros(G.cells.num,1);

G = computeVEM2DGeometry(G,f,1,1);

sourceCoords = [.3,.2];
source = sum(bsxfun(@minus, G.cells.centroids, sourceCoords).^2,2);
source = find(source == min(source));
source = source(1);

sinkCoords = [.5,.5];
sink = sum(bsxfun(@minus, G.cells.centroids, sinkCoords).^2,2);
sink = find(sink == min(sink));
sink = sink(1);

f = @(X) zeros(size(X,1),1);

Q = 10;
src = [];
% src = addSource(src, source, Q);
src = addSource(src, sink, -Q);

tol = 1e-6;
boundaryEdges = find(G.faces.neighbors(:,1) == 0 | G.faces.neighbors(:,2) == 0);
isExternal = abs(G.faces.centroids(boundaryEdges,1)) < tol | ...
             abs(G.faces.centroids(boundaryEdges,1) - gridLim(1)) < tol | ...
             abs(G.faces.centroids(boundaryEdges,2)) < tol |...
             abs(G.faces.centroids(boundaryEdges,2) - gridLim(2)) < tol;
isInternal = ~isExternal;

bc = VEM_addBC(G, [], boundaryEdges(isExternal), 'pressure', 0);
% bc = VEM_addBC(G, bc, boundaryEdges(isInternal), 'flux', 0);
            
sol1 = VEM2D_v3(G,f,1,bc,'src', src);
sol1 = cellAverages(G,sol1);

sol2 = VEM2D_v3(G,f,2,bc, 'src', src);

figure;
plotCellData(G, sol1.cellAverages);
plotFaces(G,boundaryEdges(isInternal),'EdgeColor', 'r')
axis([0,1,0,1])
colorbar;
figure;
plotCellData(G, sol2.cellMoments)
plotFaces(G,boundaryEdges(isInternal),'EdgeColor', 'r')
axis([0,1,0,1])
colorbar;
