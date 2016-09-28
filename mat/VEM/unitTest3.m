require mpfa

clear; close all;
tol = 1e-6;

%%  2D 1ST AND 2ND ORDER

G = computeVEMGeometry(unitSquare([10,10],[1,1]));


rot = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];
K = diag(100*milli*darcy*rand(2,1),0);
R = rot(rand(1)*2*pi);
K = R'*K*R;
rock.perm = repmat(K([1,2,4]), [G.cells.num,1]);

fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);

p = 1e4;

bf = boundaryFaces(G);

e = abs(G.faces.centroids(bf,1)-1) < tol;
w = abs(G.faces.centroids(bf,1)) < tol;
bc = addBC([], bf(w), 'pressure', p);
bc = addBC(bc, bf(e), 'pressure', 0);
bc = addBC(bc, bf(~e & ~w), 'flux', 0);

T = computeMultiPointTrans(G, rock);
stateMPFA = incompMPFA(state, G, T, fluid, 'bc', bc);

tic

S = computeVirtualIP(G, rock, 1);
stateVEM = incompVEM(state, G, S, fluid, 'bc', bc);

pErr = norm(stateVEM.pressure - stateMPFA.pressure)/norm(stateMPFA.pressure);
fErr = norm(stateVEM.flux - stateMPFA.flux)/norm(stateMPFA.flux);
fprintf('2D 1st order: \t pressure error: \t %5.2f %%\t', pErr*100);
fprintf('flux error: \t %5.2f%%\t ', fErr*100);

toc

tic

S = computeVirtualIP(G, rock, 2);
stateVEM = incompVEM    (state, G, S, fluid, 'bc', bc);

pErr = norm(stateVEM.pressure - stateMPFA.pressure)/norm(stateMPFA.pressure);
fErr = norm(stateVEM.flux - stateMPFA.flux)/norm(stateMPFA.flux);
fprintf('2D 2nd order: \t pressure error: \t %5.2f %%\t', pErr*100);
fprintf('flux error: \t %5.2f%%\t ', fErr*100);

toc

%%  3D 1ST AND 2ND ORDER

G = computeVEMGeometry(voronoiCube(100,[1,1,1]));

rot = @(theta, eta, psi) [1 0 0; 0 cos(theta), -sin(theta); 0 sin(theta) cos(theta)]...
                        *[cos(eta) 0 sin(eta); 0 1 0; -sin(eta) 0 cos(eta)         ]...
                        *[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1         ];

K = diag(100*milli*darcy*rand(1,3),0);
R = rot(rand(1)*2*pi, rand(1)*2*pi, rand(1)*2*pi);
K = R'*K*R;
rock.perm = repmat(K([1,2,3,5,6,9]), [G.cells.num, 1]);

% rock.perm = ones(G.cells.num,1);

fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);
state = initState(G, [], 0);

p = 1e4;

bf = boundaryFaces(G);

e = abs(G.faces.centroids(bf,1)-1) < tol;
w = abs(G.faces.centroids(bf,1)) < tol;
bc = addBC([], bf(w), 'pressure', p);
bc = addBC(bc, bf(e), 'pressure', 0);
bc = addBC(bc, bf(~e & ~w), 'flux', 0);

T = computeMultiPointTrans(G, rock);
stateMPFA = incompMPFA(state, G, T, fluid, 'bc', bc);

tic

S = computeVirtualIP(G, rock, 1);
stateVEM = incompVEM(state, G, S, fluid, 'bc', bc, 'facePressure', true);

pErr = norm(stateVEM.pressure - stateMPFA.pressure)/norm(stateMPFA.pressure);
fErr = norm(stateVEM.flux - stateMPFA.flux)/norm(stateMPFA.flux);
fprintf('3D 1st order: \t pressure error: \t %5.2f %%\t', pErr*100);
fprintf('flux error: \t %5.2f%%\t ', fErr*100);

toc

tic

S = computeVirtualIP(G, rock, 2);
stateVEM = incompVEM(state, G, S, fluid, 'bc', bc);

pErr = norm(stateVEM.pressure - stateMPFA.pressure)/norm(stateMPFA.pressure);
fErr = norm(stateVEM.flux - stateMPFA.flux)/norm(stateMPFA.flux);
fprintf('3D 2nd order: \t pressure error: \t %5.2f %%\t', pErr*100);
fprintf('flux error: \t %5.2f%%\t ', fErr*100);

toc