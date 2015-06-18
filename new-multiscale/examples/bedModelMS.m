pth = getDatasetPath('bedmodel2');

grdecl = readGRDECL(fullfile(pth, 'BedModel2.grdecl'));
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

G = processGRDECL(grdecl);
G = computeGeometry(G);

rock = grdecl2Rock(grdecl, G.cells.indexMap);

%%
mrstModule add mrst-gui
figure;
plotToolbar(G, rock);
axis tight off
view(30, 30);

figure;
plotToolbar(G, grdecl.SATNUM(G.cells.indexMap));
axis tight off
view(30, 30);

%%
xf = G.faces.centroids(:, 1);
left = find(abs(xf - min(xf)) < 1e-4);
right = find(abs(xf - max(xf)) < 1e-4);

bc = [];
bc = addBC(bc, left, 'pressure', 1);
bc = addBC(bc, right, 'pressure', 0);


figure;
plotGrid(G, 'facecolor', 'none');
plotFaces(G, left, 'facecolor', 'r')
plotFaces(G, right, 'facecolor', 'b')
view(30, 30)
%%
fluid = initSingleFluid('rho', 1, 'mu', 1);
T = computeTrans(G, rock);
state = initResSol(G, 0);

%%
ref = incompTPFA(state, G, T, fluid, 'bc', bc);
%% Grid of roughly same dimensions
dims = max(G.nodes.coords);
Gc = cartGrid([30, 30, 60], dims);
Gc = computeGeometry(Gc);
%% Make coarse grid
p = partitionUI(Gc, [6, 6, 5]);
p = reshape(p, Gc.cartDims);
p = sampleFromBox(G, p);
% p = partitionMETIS(G, T, 200);
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegion(CG, 'edgeBoundaryCenters', true, 'adjustCenters', true);
%%
A = getIncomp1PhMatrix(G, T, state, fluid);
basis = getMultiscaleBasis(CG, A, 'type', 'msrsb', 'iterations', 150);

%%
ms = incompMultiscale(state, CG, T, fluid, basis, 'bc', bc);

%%
figure;
plotToolbar(G, ms)
title('MS')

figure;
plotToolbar(G, ref)
title('TPFA')

%%
for i = 1:CG.cells.num
    basis.B(CG.cells.centers(i), i)
end