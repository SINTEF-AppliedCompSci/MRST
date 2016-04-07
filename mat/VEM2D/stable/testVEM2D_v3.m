clc; clear all; close all;

addpath('../')
addpath('../../');
addpath('~/Documents/master/pebiGridding/voronoi2D/')
run('../../startup.m');

n = 10;

% G = cartGrid([n,n], [1,1]);
G = unitSquare(n,n);
f = @(X) 0*ones(size(X,1),1);
g = @(X) X(:,1);
% f = zeros(G.cells.num,1);
G = computeVEM2DGeometry(G,f,1,1);

boundaryEdges = find(G.faces.neighbors(:,1) == 0 | G.faces.neighbors(:,2) == 0);

bc = VEM_addBC(G, [], boundaryEdges, 'pressure', g);
            
sol1 = VEM2D_v3(G,f,1,bc);
sol1 = cellAverages(G,sol1);

sol2 = VEM2D_v3(G,f,2,bc);

figure;
plotCellData(G, sol1.cellAverages);
colorbar;
figure;
plotCellData(G, sol2.cellMoments)
colorbar;

u = [g(G.nodes.coords); g(G.faces.centroids); polygonInt_v2(G,1:G.cells.num,g,7)];
U1 = sol1.nodeValues;
U2 = [sol2.nodeValues; sol2.edgeValues; sol2.cellMoments];

norm(u(1:G.nodes.num) - U1)
norm(u-U2)

norm(sol1.nodeValues-sol2.nodeValues)