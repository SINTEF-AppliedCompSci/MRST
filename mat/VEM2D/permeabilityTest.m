clc; clear; close all;
mrstModule add mpfa
mrstModule add mimetic

gravity on

%%  Create gird, set random permeability.

n = 17;
G = cartGrid([n,n],[1,1]);
G = sortEdges(G);
G = computeVEM2DGeometry(G);

c = 3;

switch c
    case 1
        rock.perm = repmat([1,1],         G.cells.num, 1);
    case 2
        rock.perm = repmat([1,1]*1e-12,      G.cells.num, 1);
    case 3
        rock.perm = repmat([10,1]*1e-12, G.cells.num, 1);
    case 4
        rock.perm = repmat([1,10]*1e-12, G.cells.num, 1);
    case 5
        p = gaussianField(G.cartDims, [0.2 0.4], [11 3], 2.5);
        K = p.^3.*(1e-5)^2./(0.81*72*(1-p).^2);
        rock.perm = K(:);
end

fluid = initSingleFluid('mu', 100, 'rho', 1);

%% Set boundary conditions.

h = .7;
bf = boundaryFaces(G);
bfw = G.faces.centroids(bf, 1) == 0 & G.faces.centroids(bf,2) < h;
bfe = G.faces.centroids(bf, 1) == 1 & G.faces.centroids(bf,2) > 1-h;

bcVEM = [];
bcVEM = VEM2D_addBC(bcVEM, G, bf(bfw), 'pressure', 1e3);
bcVEM = VEM2D_addBC(bcVEM, G, bf(bfe), 'pressure', .0);
bcVEM = VEM2D_addBC(bcVEM, G, bf(~bfw & ~bfe), 'flux', 0);

bcMRST = [];
bcMRST = addBC(bcMRST, bf(bfw), 'pressure', 1e3);
bcMRST = addBC(bcMRST, bf(bfe), 'pressure', 0);
bcMRST = addBC(bcMRST, bf(~bfw & ~bfe), 'flux', 0);

%% Initialize states, compute transmissibilities

state = initState(G, [], 0);
T = computeTrans(G, rock);
TMPFA = computeMultiPointTrans(G, rock);
S = computeMimeticIP(G, rock);

%%  Solve

solVEM = VEM2D(state, G, rock, fluid, 2, 'bc', bcVEM, 'cellPressure', true);
solTPFA = incompTPFA(state, G, T, fluid, 'bc', bcMRST);
solMPFA = incompMPFA(state, G, TMPFA, fluid, 'bc', bcMRST);
solMFD  = solveIncompFlow(state, G, S, fluid, 'bc', bcMRST);

%% Plot results

figure;
subplot(2,2,1)
plotCellData(G,solVEM.pressure)
colorbar
axis equal
axis([0,1,0,1])

subplot(2,2,2)
plotCellData(G, solTPFA.pressure);
colorbar
axis equal
axis([0,1,0,1])

subplot(2,2,3)
plotCellData(G, solMPFA.pressure);
colorbar
axis equal
axis([0,1,0,1])

subplot(2,2,4)
plotCellData(G, solMFD.pressure);
colorbar
axis equal
axis([0,1,0,1])

colormap(parula)

%% Calculate and plot error

errTPFA = abs(solTPFA.pressure-solVEM.pressure);
errMPFA = abs(solMPFA.pressure-solVEM.pressure);
errMFD = abs(solMFD.pressure-solVEM.pressure);

figure;
subplot(1,3,1)
plotCellData(G,errTPFA)
colorbar
axis equal
axis([0,1,0,1])

subplot(1,3,2)
plotCellData(G, errMPFA);
colorbar
axis equal
axis([0,1,0,1])

subplot(1,3,3)
plotCellData(G, errMFD);
colorbar
axis equal
axis([0,1,0,1])

colormap(parula)

errTPFANorm = norm(errTPFA,2)/max(abs(solTPFA.pressure));
errMPFANorm = norm(errMPFA,2)/max(abs(solMPFA.pressure));
errMFDNorm = norm(errMFD,2)/max(abs(solMFD.pressure));

fprintf('Relative difference: \n-----------------\n');
fprintf('TPFA \t %5.2f %%\n', errTPFANorm*100);
fprintf('MPFA \t %5.2f %%\n', errMPFANorm*100);
fprintf('MFD  \t %5.2f %%\n', errMFDNorm*100);