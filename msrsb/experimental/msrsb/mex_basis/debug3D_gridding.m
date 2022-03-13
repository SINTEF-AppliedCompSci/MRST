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

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
