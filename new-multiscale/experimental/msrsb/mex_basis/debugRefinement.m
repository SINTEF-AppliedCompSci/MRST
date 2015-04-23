mrstModule add multiscale-devel coarsegrid 

if 0
    G = cartGrid([2 1]);
    p = [1; 2];
elseif 0
    G = cartGrid([2 2]);
    p = partitionUI(G, [2 2]);
elseif 0
    G = cartGrid([30 90]);
    p = partitionUI(G, [3 3]);
    refined = p == 5;
    p(refined) = (1:sum(refined))' + 9;
else
    G = cartGrid([100 200]);
    cdims = ceil(G.cartDims./[10 20]);
    p = partitionUI(G, cdims);
%     p = partitionUI(G, [3 3]);
    
    [ii, jj] = gridLogicalIndices(G);
    refined = ii < 6 & jj < 6;
    p(refined) = (1:sum(refined))' + 9;
end
close all
plotToolbar(G, mod(p, 13))

G = computeGeometry(G);

rock.perm = repmat(500*milli*darcy, G.cells.num, 1);
A = getIncomp1PhMatrix(G, computeTrans(G, rock));




p = processPartition(G, p);
p = compressPartition(p);

CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegion(CG);

%%
sg = setupGridsForMex(A, CG);
active = true(G.cells.num, 1);

% active = refined;


msdir = mrstPath('query', 'multiscale-devel');
delete(fullfile(msdir, 'mex_basis', 'mex_iteratedJacobiBasisFaster.mex*'))
rehash

[I, I_comp] = mbasisSmoothed(sg, 'tolerance',2, 'maxiter', 1000, 'omega', .66, 'active', active);

clear tmp
tmp.sum = sum(I, 2);
tmp.I = I;

close all
plotToolbar(G, tmp)
colorbar
%%
% sg = setupGridsForMex(A, CG);

close all
for i = 1:CG.cells.num
    figure(1); clf;
    
    c = sg.cells{i} + 1;
    cbnd = c(sg.isBnd{i} == 1);
    
    c = CG.cells.interaction{i};
    
    plotGrid(G, c, 'FaceColor', 'red');
%     plotGrid(G, setdiff(c, cbnd), 'FaceColor', 'red');
    plotGrid(G, cbnd, 'FaceColor', 'blue', 'FaceAlpha', .2);
    outlineCoarseGrid(G, p);
    pause()
end
