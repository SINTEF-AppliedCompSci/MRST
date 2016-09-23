clc; clear all; close all;

%%  2D 1ST ORDER

fprintf('Test 1 of 4, 2D 1st order ... ')

tic

G = computeVEMGeometry(unitSquare([10,10],[1,1]));

gD = @(x) x(:,1) + 10*x(:,2);

bf = boundaryFaces(G);
bc = addBC([], bf, 'pressure', gD(G.faces.centroids(bf,:)));

rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);

S = computeVirtualIP(G, rock, 1);
state = incompVEM(state, G, S, fluid, 'bc', bc);

err = norm(state.nodePressure - gD(G.nodes.coords));
fprintf('Error: \t %.2d.\t', err);

toc;

%%  2D 2ND ORDER

fprintf('Test 2 of 4, 2D 2nd order ... ')

tic

G = computeVEMGeometry(unitSquare([10,10],[1,1]));
G = sortEdges(G);
gD = @(x) x(:,1).^2 - x(:,2).^2 + x(:,1).*x(:,2);
bf = boundaryFaces(G);
bc = addBCFunc([], bf, 'pressure', gD);
rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);
S = computeVirtualIP(G, rock, 2);
state = incompVEM(state, G, S, fluid, 'bc', bc);

p = [gD(G.nodes.coords); gD(G.faces.centroids); ...
     polygonInt(G, 1:G.cells.num, gD, 2)./G.cells.volumes];
P = [state.nodePressure; state.facePressure; state.cellPressure];

err = norm(p-P);
fprintf('Error: \t %.2d.\t', err);

toc;

%%  3D 1ST ORDER

fprintf('Test 3 of 4, 3D 1st order ... ')

tic

G = computeVEMGeometry(voronoiCube(100,[1,1,1]));

gD = @(x) x(:,1) + x(:,2) + x(:,3);

bf = boundaryFaces(G);
bc = addBCFunc([], bf, 'pressure', gD);

rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);

S = computeVirtualIP(G, rock, 1);
state = incompVEM(state, G, S, fluid, 'bc', bc);

err = norm(state.nodePressure - gD(G.nodes.coords));
fprintf('Error: \t %.2d.\t', err);

toc;

%%  3D 2ND ORDER

fprintf('Test 4 of 4, 3D 2nd order ... ')

tic

G = computeVEMGeometry(voronoiCube(100,[1,1,1]));

gD = @(x) x(:,1).^2 - x(:,2).^2;

bf = boundaryFaces(G);
bc = addBCFunc([], bf, 'pressure', gD);

rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);

S = computeVirtualIP(G, rock, 2);
state = incompVEM(state, G, S, fluid, 'bc', bc);

p = [gD(G.nodes.coords); gD(G.edges.centroids); ...
     polygonInt3D(G, 1:G.faces.num, gD, 2)./G.faces.areas; ...
     polyhedronInt(G, 1:G.cells.num, gD, 2)./G.cells.volumes];
P = [state.nodePressure; state.edgePressure; state.facePressure; state.cellPressure];

err = norm(p-P);
fprintf('Error: \t %.2d.\t', err);

toc;
