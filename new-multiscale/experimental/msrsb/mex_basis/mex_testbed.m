mrstModule add mex spe10 multiscale-devel coarsegrid mrst-gui libgeometry
if 0
    [G, W, rock] = SPE10_setup(1);
    iscart = true;
elseif 0
    load /data/norne
    iscart = false;
else
    G = cartGrid([100, 100, 10], [3 3 2]);
    rock = struct('perm', ones(G.cells.num, 1));
    % rock = SPE10_rock(1:G.cartDims(1), 1:G.cartDims(2), 1:G.cartDims(3));
    % rock.perm = convertFrom(rock.perm(:, 1), milli*darcy);
    G = mcomputeGeometry(G);
    iscart = false;
end
fluid = initSingleFluid('mu', 1, 'rho', 1);

state = initResSol(G, 0);
T = computeTrans(G, rock);
% T(T<1e-12) = 0;
% T(52582) = 0;

state = incompTPFA(state, G, T, fluid, 'MatrixOutput', true, 'LinSolve', @(A, x) 0*x);

cdims = ceil(G.cartDims./[10, 10, 10]);
% cdims = ceil(G.cartDims./[5 5 5]);
% cdims = ceil(G.cartDims./[15 15 1]);

if 0
    p = partitionMETIS(G, T, prod(cdims));
    iscart = false;
else
    p = partitionUI(G, cdims);
end
p = processPartition(G, p);
p = compressPartition(p);
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);

if iscart
    CG = storeInteractionRegionCart(CG);
else
    CG = storeInteractionRegion(CG);
end
A = state.A;
%%
clc
msdir = mrstPath('query', 'multiscale-devel');
tic()
sg = setupGridsForMex(A, CG);
toc()
%%
delete(fullfile(msdir, 'mex_basis', 'mex_iteratedJacobiBasisFaster.mex*'))

rehash
I = mbasisSmoothed(sg, 'tolerance', 2, 'maxiter', 1, 'omega', .66);


active = ismember(CG.partition, 5);
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
clc
msdir = mrstPath('query', 'multiscale-devel');

delete(fullfile(msdir, 'mex_basis', 'mex_iteratedJacobiBasis.mex*'))
rehash
[I, ii, overlap] = miteratedJacobiBasis(A, CG);
figure
clear tmp;
tmp.I = I;
tmp.sumi = sum(I, 2);
plotToolbar(G, tmp)
outlineCoarseGrid(G, CG.partition)
%%
b = rand(G.cells.num, 1);

x = simpleMultiscaleGMRES(A, b, I);
%%
tic
I = iteratedJacobiBasis(A, CG, 'iterations', 100, 'autostop', false, 'verbose', true);
toc
close all
clear tmp;
tmp.I = I;
tmp.sumi = sum(I, 2);
plotToolbar(G, tmp)
 outlineCoarseGrid(G, CG.partition)



%%
tmp = nan(G.cells.num, CG.cells.num);
for i = 1:CG.cells.num
    tmp(CG.cells.interaction{i}, i) = 1;
end
close all
plotToolbar(G, tmp);
%%
tmp = nan(G.cells.num, CG.cells.num+1);
for i = 1:CG.cells.num
    tmp(sg.cells{i}(sg.isBnd{i}) + 1, i) = 1;
end
tmp(:, end) = sum(~isnan(tmp), 2);
figure;
plotToolbar(G, tmp);
%%
bad = false(G.cells.num, 1);
for i = 1:G.cells.num
    t = sum(abs(A(i, :)) > 0);
    if t < 2
        bad(i) = true;
    end
end
