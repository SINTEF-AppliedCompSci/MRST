clc; clear all; close all;

addpath('../VEM3D/');
mrstModule add mimetic

n = 100;
k = 2;
G = cartGrid([n,n]);


G = computeVEM2DGeometry(G);
G = sortEdges(G);

state = initState(G, [], 0);
rock.perm = ones(G.cells.num,1);
fluid = initSingleFluid('mu', 1, 'rho', 1);
src = addSource([], [1,2], 1);

tic;
S = computeVirtualIP(G, rock, fluid, k);
toc

tic;
SMim = computeMimeticIP(G, rock);
toc


state = incompVEM(state, G, S, fluid, 'src', src);