clc; clear all; close all;

addpath('../VEM3D/');

G = cartGrid([300,300]);

G = computeVEM2DGeometry(G);
G = sortEdges(G);

state = initState(G, [], 0);
rock.perm = ones(G.cells.num,1);
fluid = initSingleFluid('mu', 1, 'rho', 1);

tic;
A = computeVirtualIP(G, rock, fluid, 1);
toc

tic;
S = computeMimeticIP(G, rock);
toc
% 
% tic;
% state = incompVEM2D(state, G, rock, fluid, 1);
% toc