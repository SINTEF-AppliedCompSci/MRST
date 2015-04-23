mrstModule add mex spe10 multiscale-devel coarsegrid mrst-gui libgeometry

G = cartGrid([30, 30, 30], [10 10 1]);
rock = struct('perm', ones(G.cells.num, 1));
G = mcomputeGeometry(G);
iscart = true;

fluid = initSingleFluid('mu', 1, 'rho', 1);

state = initResSol(G, 0);
T = computeTrans(G, rock);

cdims = ceil(G.cartDims./[10, 10, 10]);

p = partitionUI(G, cdims);

p = processPartition(G, p);
p = compressPartition(p);
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);


% CG = storeInteractionRegionCart(CG);
CG = storeInteractionRegion(CG);

A = getIncomp1PhMatrix(G, T);

clc
msdir = mrstPath('query', 'multiscale-devel');
tic()
sg = setupGridsForMex(A, CG);
toc()

figure; plotGrid(G, sg.globalBoundary+1)
%%
delete(fullfile(msdir, 'mex_basis', 'mex_iteratedJacobiBasisFaster.mex*'))

rehash
I = mbasisSmoothed(sg, 'tolerance', 2, 'maxiter', 1, 'omega', .66);


% active = ismember(CG.partition, 5);
active = true(G.cells.num, 1);
I = mbasisSmoothed(sg, 'tolerance', 2, 'maxiter', 1000, 'omega', .66, 'active', active);

% [I, ii, overlap] = miteratedJacobiBasisFaster(A, CG);
figure
clear tmp;
tmp.I = I;
tmp.sumi = sum(I, 2);

[min(tmp.sumi), max(tmp.sumi)]
plotToolbar(G, tmp)
axis tight off

% outlineCoarseGrid(G, CG.partition)
colorbar


%%
tmp = nan(G.cells.num, CG.cells.num);
for i = 1:CG.cells.num
    tmp(CG.cells.interaction{i}, i) = 1;
end
close all
plotToolbar(G, tmp, 'edgec', 'k');
