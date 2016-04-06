clc; clear all; close all;

addpath('../../');
addpath('~/Documents/master/pebiGridding/voronoi2D/')
run('../../startup.m');

n = 10;

% G = cartGrid([n,n], [1,1]);
G = unitSquare(n,n);
% f = @(X) -2*ones(size(X,1),1);
f = @(X) X(:,1) + 100*X(:,2);
% f = zeros(G.cells.num,1);
G = computeVEM2DGeometry(G,f,1,1);

boundaryEdges = find(G.faces.neighbors(:,1) == 0 | G.faces.neighbors(:,2) == 0);

bc = addBC([], boundaryEdges, 'pressure', 0);
            
sol1 = VEM2D_v3(G,f,1,bc);
sol1 = cellAverages(G,sol1);

sol2 = VEM2D_v3(G,f,2,bc);

figure;
plotCellData(G, sol1.cellAverages);
colorbar;
figure;
plotCellData(G, sol2.cellMoments)
colorbar;

u = 

norm(sol1.nodeValues-sol2.nodeValues)