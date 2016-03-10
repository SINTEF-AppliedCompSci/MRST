%   TEST 2: unit suqare 3D
clc; clear; close all;

addpath('../'); addpath('../VEM3D/');

%   TEST 4: Finite difference 3D

nx = 1; ny = 1; nz = 1;
G = cartGrid([nx,ny,nz],[2,2,2]);
X = G.nodes.coords;
X = bsxfun(@minus, X, [1,1,1]);
G.nodes.coords = X;

k = 1;
f = @(X) zeros(size(X,1),1);
G = computeVEMGeometry(G,f,k);

l1 = sqrt(8/9);
l4 = sqrt(8/27);

a1 = sqrt(16/3);
a4 = sqrt(8/3);

h = G.cells.diameters(1);
      
AK_FEM = retrieveFEMMatrix(3,1);


alp = 0:0.1:1;

alpha = repmat([6, 6, 6, 9],G.cells.num,1);

boundaryEdges = find((G.faces.neighbors(:,1) == 0) + (G.faces.neighbors(:,2) == 0));
gD = @(X) zeros(size(X,1),1);
bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryEdges}}, 'bcType', {{'dir'}});

[~, A_VEM,~,G] = VEM3D(G,f,bc,1,alpha);

AK_VEM = G.cells.AK;
err = norm(AK_FEM - AK_VEM, 'fro');
fprintf('Result from unit square comparison: \t  %d \n\n', err);