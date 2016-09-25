clc; clear; close all;
tol = 1e-6;

%%  2D 1ST ORDER

fprintf('Test 1 of 4, 2D 1st order ... ')

tic

G = computeVEMGeometry(unitSquare([10,10],[1,1]));

rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);


Q = 1;
bf = boundaryFaces(G);
e  = abs(G.faces.centroids(bf,1)-1) < tol;
bc = addBC([], bf(~e), 'pressure', 0);
bc = addBC(bc, bf( e), 'flux'    , -Q*G.faces.areas(bf(e))/(sum(G.faces.areas(bf(e)))));

C = [.25,.5];
d = sum(bsxfun(@minus, G.cells.centroids, C).^2,2);
srcCell = find(d == min(d));
src = addSource([], srcCell(1), Q);

S = computeVirtualIP(G, rock, 1);
stateVEM = incompVEM    (state, G, S, fluid, 'bc', bc, 'src', src);
stateVEM = calculateCellPressure(stateVEM, G, S);

S = computeMimeticIP(G, rock);
stateMFD = incompMimetic(state, G, S, fluid, 'bc', bc, 'src', src);


err = norm(stateVEM.cellPressure - stateMFD.pressure)/norm(stateMFD.pressure);
fprintf('Error: \t %.2f %%\t', err*100);

subplot(1,2,1);
plotCellData(G, stateVEM.cellPressure);
colorbar;

subplot(1,2,2);
plotCellData(G, stateMFD.pressure);
colorbar;

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
bc = addBCVEM(bc, bf( e), 'flux'    , gN);

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
rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);


gD = @(x) x(:,1) + x(:,2) + x(:,3);
bf = boundaryFaces(G);
bc = addBCVEM([], bf, 'pressure', gD);

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
bc = addBCVEM([], bf, 'pressure', gD);

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
