clc; clear all; close all;
addpath('../');
addpath('../VEM2D/');

%   TEST 1: unit square

G = cartGrid([2,2], [2,2]);

X = G.nodes.coords;
X = bsxfun(@minus, X, [1,1]);
G.nodes.coords = X;

f = @(X) X(:,1);

k = 1;

G = computeVEM2DGeometry(G,f,k,2/3);

load('FEM2D_1st_sq');
AK_VEM = G.cells.AK{1};

err = norm(AK_FEM - AK_VEM, 'fro');

fprintf('Result from unit square comparison: \t  %d \n\n', err);

k = 2;
G = computeVEM2DGeometry(G,f,k,0);
AK_VEM = G.cells.AK{1};
SK = G.cells.SK{1};

load('FEM2D_2nd_sq');

err = norm(AK_FEM - AK_VEM, 'fro');

fprintf('Result from unit square comparison: \t  %d \n\n', err);