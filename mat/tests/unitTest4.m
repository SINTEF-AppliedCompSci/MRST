require mimetic

clear; close all;
tol = 1e-6;

%%  2D 1ST AND 2ND ORDER

G = computeVEMGeometry(cartGrid([11,11,11], [1,1,1]));

rock.perm = 100*milli*darcy*ones([G.cells.num,1]);

fluid      = initSingleFluid('mu' , 1*centi*poise, ...
                             'rho', 1000*kilogram/meter^3);

bc = addBC([], boundaryFaces(G), 'pressure', 0);

C = [.5,.5];
d = sum(bsxfun(@minus, G.cells.centroids(:,[1,2]), C).^2,2);
wc = find(d == min(d));

W = addWell([], G, rock, wc, 'val', 1e4, 'radius', .0001, 'type', 'bhp');
state = initState(G, W, 0);

tic

S = computeVirtualIP(G, rock, 2);
stateVEM = incompVEM(state, G, S, fluid, 'bc', bc, 'wells', W);

toc

S = computeMimeticIP(G, rock);
stateMFD = incompMimetic(state, G, S, fluid, 'bc', bc, 'wells', W);

T = computeMultiPointTrans(G, rock);
stateMPFA = incompMPFA(state, G, T, fluid, 'bc', bc, 'wells', W);

dd = norm(stateVEM.pressure - stateMFD.pressure)/norm(stateMFD.pressure);
fprintf('Difference MFD: \t %.2f %%\n', dd*100);

dd = norm(stateVEM.pressure - stateMPFA.pressure)/norm(stateMPFA.pressure);
fprintf('Difference MPFA: \t %.2f %%\n', dd*100);

%%

c = G.cells.centroids(:,2) < .5 | (d == min(d));

subplot(1,3,1)
plotCellData(G, stateVEM.pressure(c), c);
set(gca, 'zdir', 'normal')
axis equal off;
view([200,4]);
colorbar

subplot(1,3,2)
plotCellData(G, stateMFD.pressure(c), c);
set(gca, 'zdir', 'normal')
axis equal off;
view([200,4]);
colorbar

subplot(1,3,3)
plotCellData(G, abs(stateMFD.pressure(c)-stateVEM.pressure(c)), c);
set(gca, 'zdir', 'normal')
axis equal off;
view([200,4]);
colorbar

%%