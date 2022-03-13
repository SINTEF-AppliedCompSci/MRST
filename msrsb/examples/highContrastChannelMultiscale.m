%% Problem comparing different multiscale methods on an idealized channel
% Extensive numerical experiments have shown that contemporary multiscale
% methods provide approximate solutions of good quality for highly
% heterogeneous media. However, cases with large permeability contrasts are
% generally challenging and it is not difficult to construct pathological
% test cases on which a particular method fails to produce accurate
% solutions. Diagonal channel. Kippe et al. proposed a simple and
% illuminating example consisting of a narrow high-permeable channel in a
% low-permeable background, where the channel is aligned with the diagonal
% direction. 
%
% This example corresponds to Example 4.3.3 in
%  A multiscale restriction-smoothed basis method for high contrast porous
%  media represented on unstructured grids. J. Comput. Phys, Vol. 304, pp.
%  46-71, 2016. DOI: 10.1016/j.jcp.2015.10.010

% Load modules
mrstModule add msrsb incomp coarsegrid mrst-gui msfvm matlab_bgl

%% We set up the domain
% The domain consists of 65 by 65 fine cells and a permeability field that
% is equal to 10 mD everywhere except for a narrow channel that crosses the
% domain diagonally, where the permeability is 1000 mD.

dims = [65, 65];
pdims = [1000 1000]*meter;
G = cartGrid(dims, pdims);
G = computeGeometry(G);

% Alternate intersection can be inserted. This channel crosses the blocks
% diagonally and is not a problem for either finite-volume multiscale
% scheme.
% x1 = repmat([0, 0, 0], G.cells.num, 1);
% x2 = repmat([1, 1, 0], G.cells.num, 1);
x1 = repmat([0.0, 0.05, 0], G.cells.num, 1);
x2 = repmat([0.95, 1.0, 0], G.cells.num, 1);

% Insert channel based on distance
x0 = [bsxfun(@rdivide, G.cells.centroids, pdims), zeros(G.cells.num, 1)];
dist = sqrt(sum(cross(x2 - x1, x1 - x0).^2, 2)./sum((x2 - x1).^2, 2));
isChannel = dist < 0.02;
% Set up rock structure
perm = repmat(10*milli*darcy, G.cells.num, 1);
perm(isChannel) = 1*darcy;
rock = makeRock(G, perm, 1);
% Setup unit fluid model
fluid = initSimpleFluid('rho', [1, 1], 'mu', [1, 1], 'n', [1 1]);
% Set up a coarse grid
coarsedims = [8, 8];
p = partitionUI(G, coarsedims);
% Plot the permeability field and coares grid
figure(1); clf
plotToolbar(G, rock.perm/darcy());
colorbar
axis equal tight off
outlineCoarseGrid(G, p, 'w')

%% Solve fine-scale problem
% We simulate 1/4 pore volume injected, corresponding to the channel being
% completely filled. Some fluid will infiltrate into the low permeability
% region, which is where we will see the approximation error made by the
% multiscale methods clearly.
time = 1*year;
Nt = 20;
pv = poreVolume(G, rock);
% Set rate and wells in opposite corners of grid
injr = sum(pv)/(4*time);
W = [];
W = verticalWell(W, G, rock, 1, 1, [], 'type', 'rate', ...
    'val', injr, 'comp_i', [1, 0]);
W = verticalWell(W, G, rock, G.cartDims(1), G.cartDims(2), [],...
    'type', 'bhp' , 'val', 1*barsa, 'comp_i', [1, 0]);

% Solve pressure field for fine scale
state0 = initResSol(G, 0, [0, 1]);
T = computeTrans(G, rock);
ref = incompTPFA(state0, G, T, fluid, 'MatrixOutput', true, 'Wells', W);

%% Coarse grid and basis functions
% Grid structure from partition vector
CG = generateCoarseGrid(G, p);
% Add centroids / geometry information on coarse grid
CG = coarsenGeometry(CG);
% Store the support of each cell (used by multiscale basis construction)
% Set up coarse grid
blockPts = CG.cells.centroids;
% Set the interaction/support regions to be as simple as possible for MsRSB
CG = storeInteractionRegion(CG, 'adjustCenters', false, ...
                                'localTriangulation', false, ...
                                'centerOverride', blockPts);
% Set up dual grid - required for classical MsFV-method
DG = partitionUIdual(CG, coarsedims);
DG = makeExplicitDual(CG, DG);
CG.dual = DG;

% Compute basis functions, use very strict tolerance and an excessive
% amount of iterations to avoid any round-of error from basis functions not
% being converged. For this reason, we use a C-accelerated version of the
% MsRSB solver. If this fails to compile on your computer, you can change
% the 'useMex' argument to false and wait a little longer for the result.
A = getIncomp1PhMatrix(G, T);
basisfv = getMultiscaleBasis(CG, A, 'type', 'msfvm');
msfv = incompMultiscale(state0, CG, T, fluid, basisfv, 'wells', W);
basis = getMultiscaleBasis(CG, A, 'type', 'msrsb', 'useMex', true, ...
                                  'iterations', 1000,'tolerance', 1e-4);
% Create a second coarse grid that adapts to the already computed basis
% functions (see paper for more details).
basis2 = basis;
CG2 = basisToCoarseGrid(G, basis2.B);
basis2.R = controlVolumeRestriction(CG2.partition);

ms = incompMultiscale(state0, CG, T, fluid, basis, 'wells', W);
ms2 = incompMultiscale(state0, CG2, T, fluid, basis2, 'wells', W);


%% Simulate time loop
% We do not need to update the pressure for either solver, as the unit
% fluid results in no mobility changes during the simulation
dt = time/Nt;
tsolve = @(state) implicitTransport(state, G, dt, rock, fluid, 'wells', W);
for j = 1:Nt
    fprintf('%d of %d \n', j, Nt);

    disp('Solving MsRSB')
    ms = [ms; tsolve(ms(end))];
    disp('Solving MsRSB [a posteriori adaption]')
    ms2 = [ms2; tsolve(ms2(end))];
    disp('Solving MsFV')
    msfv = [msfv; tsolve(msfv(end))];
    disp('Solving TPFA')
    ref = [ref; tsolve(ref(end))];
end

%% Finally, plot the results
figure;
subplot(2, 2, 1)
plotCellData(G, ref(end).s(:, 1), 'EdgeColor', 'none')
axis equal tight off
caxis([0, 1])
title('Reference')

subplot(2, 2, 2)
plotCellData(G, msfv(end).s(:, 1), 'EdgeColor', 'none')
axis equal tight off
caxis([0, 1])
title('MsFV')
outlineCoarseGrid(G, CG.partition, 'w', 'linewidth', 1)

subplot(2, 2, 3)
plotCellData(G, ms(end).s(:, 1), 'EdgeColor', 'none')
axis equal tight off
caxis([0, 1])
title('MsRSB')
outlineCoarseGrid(G, CG.partition, 'w', 'linewidth', 1)

subplot(2, 2, 4)
plotCellData(G, ms2(end).s(:, 1), 'EdgeColor', 'none')
axis equal tight off
caxis([0, 1])
title('MsRSB (post-adapted grid)')
outlineCoarseGrid(G, CG2.partition, 'w', 'linewidth', 1)

%% Copyright notice

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
