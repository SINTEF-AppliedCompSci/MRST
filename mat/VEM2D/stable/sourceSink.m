clc; clear all; close all;

addpath('../../');
addpath('~/Documents/master/pebiGridding/voronoi2D/')
run('../../startup.m');

n = 5;

% x     = linspace(0.2, 0.8, 10);
% y     = 0.8 - 0.5*x - 0.05* sin(6*pi*x);
% fault = {[x' , y']};
% G = compositePebiGrid(1/20, [1, 1], ...
%                        'faultLines', fault, 'faultGridFactor', 1/sqrt(2));

x = linspace(.2,.8,10);
y = 1-x;
fault = {[x' , y']};
G = compositePebiGrid(1/40, [1, 1], ...
                       'faultLines', fault);
                   
                   
faultFaces = 1:G.faces.num;
faultFaces = faultFaces(G.faces.tag>0);
G = makeInternalBoundary(G, faultFaces);
f = zeros(G.cells.num,1);
G = computeVEM2DGeometry(G,f,1,1);
find(min(abs(bsxfun(@minus,G.cells.centroids,[.2,.2]))));

sourceCoords = [.4,.4];
source = sum(bsxfun(@minus, G.cells.centroids, sourceCoords).^2,2);
source = find(source == min(source));

sinkCoords = [.6,.6];
sink = sum(bsxfun(@minus, G.cells.centroids, sinkCoords).^2,2);
sink = find(sink == min(sink));

f = @(X) zeros(size(X,1),1);

Q = 10;
src = addSource([], source, Q);
src = addSource(src, sink, -Q);

boundaryEdges = find(G.faces.neighbors(:,1) == 0 | G.faces.neighbors(:,2) == 0);
isExternal = G.faces.centroids(boundaryEdges,1) == 0 | ...
        G.faces.centroids(boundaryEdges,1) == 1 | ...
        G.faces.centroids(boundaryEdges,2) == 0 | ...
        G.faces.centroids(boundaryEdges,2) == 1;
isIntnernal = ~isExternal;

bc = addBC([], boundaryEdges(isExternal), 'pressure', 0);
bc = addBC(bc, boundaryEdges(isIntnernal), 'flux', 0);
            
sol1 = VEM2D_v3(G,f,1,bc,'src', src);
sol1 = cellAverages(G,sol1);



sol2 = VEM2D_v3(G,f,2,bc, 'src', src);

figure;
plotCellData(G, sol1.cellAverages);
colorbar;
figure;
plotCellData(G, sol2.cellMoments)
colorbar;