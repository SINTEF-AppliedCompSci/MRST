clc; clear; close all;

%%  Create gird, set random permeability.
G = cartGrid([30,30],[1,1] );
G = computeVEM2DGeometry(G);
% p = gaussianField(G.cartDims, [0.2 0.4], [11 3], 2.5);
% K = p.^3.*(1e-5)^2./(0.81*72*(1-p).^2);
K = repmat(1e-10, G.cells.num, 1);
rock.perm = K;
fluid = initSingleFluid('mu', 1, 'rho', 1000);

state = initState(G, [], 0);

bf = boundaryFaces(G);

bfw = G.faces.centroids(bf, 1) == 0 & G.faces.centroids(bf,2) < .1;
bfe = G.faces.centroids(bf, 1) == 1 & G.faces.centroids(bf,2) > .9;

bc = VEM2D_addBC([], G, bf(bfw), 'pressure', 1e3);
bc = VEM2D_addBC(bc, G, bf(bfe), 'pressure', 0);
bc = VEM2D_addBC(bc, G, bf(~bfw & ~bfe), 'flux', 0);

sol = VEM2D(state, G, rock, fluid, 0, bc, 1, 'cellAverages', true);

% plotCellData(G, rock.perm);
figure;
plotVEM2D(G,sol,1)
colorbar