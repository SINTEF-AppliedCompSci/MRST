clc; clear all; close all;


G = cartGrid([5,5,5]);
G = computeVEM3DGeometry(G);

fluid = initSingleFluid('rho', 1, 'mu', 1);
rock.perm = ones(G.cells.num, 1);
state = initState(G, [], 0);
G = faceCoords(G);

state = incompVEM(state, G, rock, fluid, 1);
