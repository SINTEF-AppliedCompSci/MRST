clc; clear; close all;
tol = 1e-6;

%%  2D 1ST ORDER

fprintf('Test 1 of 4, 2D 1st order ... ')

tic

G = computeVEMGeometry(unitSquare([10,10],[1,1]));

rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
mu = 1*centi*poise;
fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);

gD = @(x) x(:,1) + 10*x(:,2);
gN = @(x) rock.perm(1)*ones(size(x,1), 1);

bf = boundaryFaces(G);
e  = abs(G.faces.centroids(bf,1)-1) < tol;
bc = addBCVEM([], bf(~e), 'pressure', gD);
bc = addBCVEM(bc, bf( e), 'flux'    , gN);

S = computeVirtualIP(G, rock, 1);
state = incompVEM(state, G, S, fluid, 'bc', bc);

errP = norm(state.nodePressure - gD(G.nodes.coords));
fprintf('Pressure error: \t %.2d. \t', errP);

toc;

%%  2D 2ND ORDER

fprintf('Test 2 of 4, 2D 2nd order ... ')

tic

G = computeVEMGeometry(unitSquare([10,10],[1,1]));

rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);

gD = @(x) x(:,1).^2 - x(:,2).^2;
gN = @(x) 2*rock.perm(1)*x(:,1);
bf = boundaryFaces(G);
e  = abs(G.faces.centroids(bf,1)-1) < tol;
bc = addBCVEM([], bf(~e), 'pressure', gD);
bc = addBCVEM(bc, bf( e), 'flux', gN);

S = computeVirtualIP(G, rock, 2);
state = incompVEM(state, G, S, fluid, 'bc', bc);

p = [gD(G.nodes.coords); gD(G.faces.centroids); ...
     polygonInt(G, 1:G.cells.num, gD, 2)./G.cells.volumes];
P = [state.nodePressure; state.facePressure; state.pressure];

errP = norm(p-P);
fprintf('Pressure error: \t %.2d. \t', errP);

toc;

%%  3D 1ST ORDER

fprintf('Test 3 of 4, 3D 1st order ... ')

tic

G = computeVEMGeometry(voronoiCube(100,[1,1,1]));
rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);


gD = @(x) x(:,1) + x(:,2) + x(:,3);
bf = boundaryFaces(G);
bc = addBCVEM([], bf, 'pressure', gD);

S = computeVirtualIP(G, rock, 1, 'trans', true);
state = incompVEM(state, G, S, fluid, 'bc', bc);

errP = norm(state.nodePressure - gD(G.nodes.coords));

fprintf('Pressure error: \t %.2d. \t', errP);

toc;

%%  3D 2ND ORDER

fprintf('Test 4 of 4, 3D 2nd order ... ')

tic

G = computeVEMGeometry(voronoiCube(100,[1,1,1]));

K = 100*milli*darcy*rand(1,3);
rock.perm = repmat(K, [G.cells.num,1]);
fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);

gD = @(x) x(:,1).^2 - K(1)/K(2)*x(:,2).^2;

bf = boundaryFaces(G);
bc = addBCVEM([], bf, 'pressure', gD);

S = computeVirtualIP(G, rock, 2, 'trans', false);
state = incompVEM(state, G, S, fluid, 'bc', bc);

p = [gD(G.nodes.coords); gD(G.edges.centroids); ...
     polygonInt3D(G, 1:G.faces.num, gD, 2)./G.faces.areas; ...
     polyhedronInt(G, 1:G.cells.num, gD, 2)./G.cells.volumes];
P = [state.nodePressure; state.edgePressure; state.facePressure; state.pressure];

errP = norm(p-P);

fprintf('Pressure error: \t %.2d. \t', errP);

toc;


