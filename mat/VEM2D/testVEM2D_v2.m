clc; clear all; close all;

addpath('../');


G = cartGrid([5,5], [1,1]);
% G = unitSquareTri([4,4],[1,1]);

f = @(X) zeros(size(X,1),1);

gD = @(X) X(:,1) + 100*X(:,2);

boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});

k = 2;
alpha = 1;

G = computeVEM2DGeometry(G,f,k, alpha);

[A,b] = VEM2D_glob_v2(G, bc, k);

U = A\b;

u = [gD([G.nodes.coords; G.faces.centroids]); polygonInt_v2(G, 1:G.cells.num, gD, 3)];

plotVEM_v2(G, U)

err = U-u;
max(abs(err))

