require mimetic

clear; close all;
tol = 1e-6;

%%  2D 1ST AND 2ND ORDER

G = computeVEMGeometry(pebiSquare([10,10],[1,1]));

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

C = [.5,.5];
d = sum(bsxfun(@minus, G.cells.centroids, C).^2,2);
srcCell = find(d == min(d));
src = addSource([], srcCell(1), Q);

S = computeMimeticIP(G, rock);
stateMFD = incompMimetic(state, G, S, fluid, 'bc', bc, 'src', src);

tr = sum(bsxfun(@minus, G.cells.centroids, C).^2,2) > 1*G.cells.diameters(srcCell(1))^2;

tic


S = computeVirtualIP(G, rock, 1, 'trans', 'mpfa');
stateVEM = incompVEM(state, G, S, fluid, 'bc', bc, 'src', src);

pErr = norm(stateVEM.pressure(tr) - stateMFD.pressure(tr))/norm(stateMFD.pressure(tr));
fErr = norm(stateVEM.flux - stateMFD.flux)/norm(stateMFD.flux);
fprintf('2D 1st order: \t pressure error: \t %5.2f %%\t', pErr*100);
fprintf('flux error: \t %5.2f%%\t ', fErr*100);

toc

tic

S = computeVirtualIP(G, rock, 2, 'trans', 'mpfa');
stateVEM = incompVEM(state, G, S, fluid, 'bc', bc, 'src', src, 'conservativeFlux', true);

pErr = norm(stateVEM.pressure(tr) - stateMFD.pressure(tr))/norm(stateMFD.pressure(tr));
fErr = norm(stateVEM.flux - stateMFD.flux)/norm(stateMFD.flux);
fprintf('2D 2nd order: \t pressure error: \t %5.2f %%\t', pErr*100);
fprintf('flux error: \t %5.2f%%\t ', fErr*100);

toc

%%  3D 1ST AND 2ND ORDER

G = computeVEMGeometry(pebiCube(200,[1,1,1]));

rot = @(theta, eta, psi) [1 0 0; 0 cos(theta), -sin(theta); 0 sin(theta) cos(theta)]...
                        *[cos(eta) 0 sin(eta); 0 1 0; -sin(eta) 0 cos(eta)         ]...
                        *[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1          ];

K1 = diag(100*milli*darcy*rand(1,3),0);
R1 = rot(rand(1)*2*pi, rand(1)*2*pi, rand(1)*2*pi);
K1 = R1'*K1*R1;
K2 = diag(100*milli*darcy*rand(1,3),0);
R2 = rot(rand(1)*2*pi, rand(1)*2*pi, rand(1)*2*pi);
K2 = R2'*K2*R2;

rock.perm = [repmat(K1([1,2,3,5,6,9]), G.cells.num/2, 1);
             repmat(K2([1,2,3,5,6,9]), G.cells.num/2, 1)];

fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);

Q = 1;
bf = boundaryFaces(G);
e  = abs(G.faces.centroids(bf,1)-1) < tol;
bc = addBC([], bf(~e), 'pressure', 0);
bc = addBC(bc, bf( e), 'flux'    , -Q*G.faces.areas(bf(e))/(sum(G.faces.areas(bf(e)))));

C = [.25, .5, .5];
d = sum(bsxfun(@minus, G.cells.centroids, C).^2,2);
srcCell = find(d == min(d));
src = addSource([], srcCell(1), Q);
C = G.cells.centroids(srcCell(1),:);

S = computeMimeticIP(G, rock);
stateMFD = incompMimetic(state, G, S, fluid, 'bc', bc, 'src', src);

tr = sum(bsxfun(@minus, G.cells.centroids, C).^2,2) > 1*G.cells.diameters(srcCell(1))^2;

tic

S = computeVirtualIP(G, rock, 1, 'trans', 'mpfa');
stateVEM = incompVEM(state, G, S, fluid, 'bc', bc, 'src', src);

pErr = norm(stateVEM.pressure(tr) - stateMFD.pressure(tr))/norm(stateMFD.pressure(tr));
fErr = norm(stateVEM.flux - stateMFD.flux)/norm(stateMFD.flux);
fprintf('3D 1st order: \t pressure error: \t %5.2f %%\t', pErr*100);
fprintf('flux error: \t %5.2f%%\t ', fErr*100);

toc

tic

S = computeVirtualIP(G, rock, 2, 'trans', 'mpfa');
stateVEM = incompVEM    (state, G, S, fluid, 'bc', bc, 'src', src);

pErr = norm(stateVEM.pressure(tr) - stateMFD.pressure(tr))/norm(stateMFD.pressure(tr));
fErr = norm(stateVEM.flux - stateMFD.flux)/norm(stateMFD.flux);
fprintf('3D 2nd order: \t pressure error: \t %5.2f %%\t', pErr*100);
fprintf('flux error: \t %5.2f%%\t ', fErr*100);

toc