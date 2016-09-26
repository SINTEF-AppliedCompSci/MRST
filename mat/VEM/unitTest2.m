require mimetic

clc; clear; close all;
tol = 1e-6;

%%  2D 1ST AND 2ND ORDER

fprintf('2D test ... \n')

G = computeVEMGeometry(unitSquare([10,10],[1,1]));

rot = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];
K    = spdiags(100*milli*darcy*rand([2*G.cells.num, 1]), 0, 2*G.cells.num, 2*G.cells.num);
R = rot(rand(1)*2*pi);
R = sparseBlockDiag(repmat(R, G.cells.num, 1), 2*ones(G.cells.num,1), 1);
K = squeezeBlockDiag(R'*K*R, 2*ones([G.cells.num,1]), 2, sum(2*G.cells.num));
K = reshape(K, 4, [])';
rock.perm = K(:, [1,2,4]);

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

S = computeMimeticIP(G, rock);
stateMFD = incompMimetic(state, G, S, fluid, 'bc', bc, 'src', src);

tic

S = computeVirtualIP(G, rock, 1);
stateVEM = incompVEM(state, G, S, fluid, 'bc', bc, 'src', src);
stateVEM = calculateCellPressure(stateVEM, G, S);

err = norm(stateVEM.cellPressure - stateMFD.pressure)/norm(stateMFD.pressure);
fprintf('1st order error: \t %5.2f %%\t', err*100);

toc

tic

S = computeVirtualIP(G, rock, 2);
stateVEM = incompVEM    (state, G, S, fluid, 'bc', bc, 'src', src);

err = norm(stateVEM.cellPressure - stateMFD.pressure)/norm(stateMFD.pressure);
fprintf('2nd order error: \t %5.2f %%\t', err*100);

toc

%%  3D 1ST AND 2ND ORDER

fprintf('3D tests... \n')

G = computeVEMGeometry(voronoiCube(500,[1,1,1]));

rot = @(theta, eta, psi) [1 0 0; 0 cos(theta), -sin(theta); 0 sin(theta) cos(theta)]...
                        *[cos(eta) 0 sin(eta); 0 1 0; -sin(eta) 0 cos(eta)         ]...
                        *[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1          ];
K = spdiags(100*milli*darcy*rand([3*G.cells.num, 1]), 0, 3*G.cells.num, 3*G.cells.num);
R = rot(rand(1)*2*pi, rand(1)*2*pi, rand(1)*2*pi);
R = sparseBlockDiag(repmat(R, G.cells.num, 1), 3*ones(G.cells.num,1), 1);
K = squeezeBlockDiag(R'*K*R, 3*ones([G.cells.num,1]), 3, sum(3*G.cells.num));
K = reshape(K, 9, [])';
rock.perm = K(:, [1,2,3,5,6,9]);
% rock.perm = 100*milli*darcy*rand(G.cells.num,3);
% rock.perm = ones(G.cells.num, 1);

fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);

Q = 1e4;
bf = boundaryFaces(G);
e  = abs(G.faces.centroids(bf,1)-1) < tol;
bc = addBC([], bf(~e), 'pressure', 0);
bc = addBC(bc, bf( e), 'flux'    , -Q*G.faces.areas(bf(e))/(sum(G.faces.areas(bf(e)))));

C = [.25, .5, .5];
d = sum(bsxfun(@minus, G.cells.centroids, C).^2,2);
srcCell = find(d == min(d));
src = addSource([], srcCell(1), Q);

S = computeMimeticIP(G, rock);
stateMFD = incompMimetic(state, G, S, fluid, 'bc', bc, 'src', src);

tr = sum(bsxfun(@minus, G.cells.centroids, C).^2,2) > 2*G.cells.diameters(srcCell(1))^2;

tic

S = computeVirtualIP(G, rock, 1);
stateVEM = incompVEM(state, G, S, fluid, 'bc', bc, 'src', src);
stateVEM = calculateCellPressure(stateVEM, G, S);

err = norm(stateVEM.cellPressure(tr) - stateMFD.pressure(tr))/norm(stateMFD.pressure(tr));
fprintf('1st order error: \t %5.2f %%\t', err*100);

toc

tic

S = computeVirtualIP(G, rock, 2);
stateVEM = incompVEM    (state, G, S, fluid, 'bc', bc, 'src', src);

err = norm(stateVEM.cellPressure(tr) - stateMFD.pressure(tr))/norm(stateMFD.pressure(tr));
fprintf('2nd order error: \t %5.2f %%\t', err*100);

toc