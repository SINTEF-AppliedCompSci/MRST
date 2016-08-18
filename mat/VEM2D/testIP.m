clc; clear all; close all;

addpath('../VEM3D/');
mrstModule add mimetic

n = 3;
k = 2;
G = cartGrid([n,n],[1,1]);

G = sortEdges(G);
G = computeVEM2DGeometry(G);

state = initState(G, [], 0);
rock.perm = ones(G.cells.num,1);
fluid = initSingleFluid('mu', 1, 'rho', 1);
src = addSource([], [1,2], 1);

f = boundaryFaces(G);

gD = @(X) X(:,2);

% bc = VEM2D_addBC([], G, f, 'pressure', gD(G.faces.centroids(f,:)));
bc = addBC([], f, 'pressure', gD(G.faces.centroids(f,:)));

tic;
S = computeVirtualIP(G, rock, fluid, k);
toc

tic;
state = incompVEM(state, G, S, fluid, 'bc', bc);
toc




plotVEM2D(G, state, 2)

% tic;
% state = VEM2D(state, G, rock, fluid, k, 'bc', bcVEM);
% toc
% 
% % plotCellData(G, state.cellPressure);