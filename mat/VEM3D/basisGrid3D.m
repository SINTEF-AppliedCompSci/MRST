clc; clear all; close all;

addpath('../VEM3D/')
% 
% G = voronoiCube(100,[1,1,1]);
% G = computeGeometry(G);
% G = mrstGridWithFullMappings(G);
% 
% intCells = find(G.cells.centroids(:,1) > .2 & ...
%                 G.cells.centroids(:,1) < .8 & ...
%                 G.cells.centroids(:,2) > .2 & ...
%                 G.cells.centroids(:,2) < .8);
% nK = numel(intCells);
% K = intCells(round(rand(1,1)*(nK-1) + 1));
% 
% nodes = G.cells.nodes(G.cells.nodePos(K):G.cells.nodePos(K+1)-1);
% P     = G.nodes.coords(nodes,:);
% 
% P(:,1) = (P(:,1) - min(P(:,1)))/(max(P(:,1))-min(P(:,1)));
% P(:,2) = (P(:,2) - min(P(:,2)))/(max(P(:,2))-min(P(:,2)));
% P(:,3) = (P(:,3) - min(P(:,3)))/(max(P(:,3))-min(P(:,3)));

P = [0   0  0; ...
     .75 0  0; ...
     1   1  0; ...
     .5  .5 1];

n = 1000;



pts = rand(n,3);

G = voronoi3D(pts, P);

plotGrid(G)