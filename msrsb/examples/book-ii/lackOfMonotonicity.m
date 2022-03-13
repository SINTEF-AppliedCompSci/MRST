%% Lack of Monotonicity (Layers of SPE10)
% This example uses a single layer from the SPE10 benchmark to illustrate
% that multiscale solutions can give nonmonotone solutions.
mrstModule add spe10 coarsegrid incomp

%%  Set up and solve the fine-scale model
% We a single-phase fluid with unit viscosity and impose a unit pressure
% drop in the y-direction.
[G, ~, rock] = getSPE10setup(85);
hT    = computeTrans(G, rock);
bc    = pside([], G, 'ymin', 1);
bc    = pside(bc, G, 'ymax', 0);
fluid = initSingleFluid('rho', 1, 'mu', 1);
state = initResSol(G, 0);

figure('Position',[300 340 780 420]);
subplot(2,2,1);
plotCellData(G, log10(rock.perm(:, 1)),'EdgeColor','none');
view(90,90), axis equal tight, title('log10(K)');

subplot(2,2,2)
state = incompTPFA(state, G, hT, fluid, 'bc', bc);
plotCellData(G, state.pressure,'EdgeColor','none');
view(90,90), axis equal tight, title('fine-scale')
colormap([1 1 1; parula(32); .7 .7 .7]);
caxis([-1/32 33/32]);
% Here, we have manipulated the colorbar and caxis so that any values out
% of the pressure range described by the boundary conditions will show up
% as white for negative values and gray for values exceeding unity.

%% Set up coarse grid
mrstModule add msrsb
cdims = ceil(G.cartDims./[5 10 5]);
p     = partitionUI(G, cdims);
CG    = generateCoarseGrid(G, p);
CG    = coarsenGeometry(CG);
CG    = storeInteractionRegionCart(CG);

%% Create dual grid and add to the coarse grid 
mrstModule add msfvm matlab_bgl
DG = partitionUIdual(CG, cdims);
DG = makeExplicitDual(CG, DG);
CG.dual = DG;

%% Compute basis functions
A        = getIncomp1PhMatrix(G, hT);
basis_sb = getMultiscaleBasis(CG, A, 'type', 'msrsb');
basis_fv = getMultiscaleBasis(CG, A, 'type', 'msfv');
bases    = {basis_sb, basis_fv};

%% Compute multiscale solutions
% Compute the MsRSB and MsFV solutions. Both schemes give multipoint
% stencils that lead to solutions containing values that are out of bounds,
% but whereas these values are barely noticable near the inflow/outflow
% boundaries for MsRSB, they are clearly evident in relatively large
% patches inside the domain for MsFV.
nb = numel(bases);
[states, reports] = deal(cell(nb, 1));
for i = 1:nb
    basis = bases{i};
    states{i} = incompMultiscale(state, CG, hT, fluid, basis, 'bc', bc);
    subplot(2,2,i+2)
    plotCellData(G, states{i}.pressure,'EdgeColor','none')
    caxis([min(state.pressure), max(state.pressure)])
    view(90,90), axis equal tight, caxis([-1/32 33/32]);
    title(bases{i}.type)
end
set(colorbar,'Position',[.925 .33 .015 .37]);

%% Solve using MS-GMRES
% As a result of the large patches with out-of-bound values, the MsFV
% method has a much higher initial residual and also converges slightly
% slower than the MsRSB method.
fn = getSmootherFunction('type', 'ilu');
for i = 1:nb
    basis = bases{i};
    [~, reports{i}] = incompMultiscale(state, CG, hT, fluid, basis, 'bc', bc,...
        'getSmoother', fn, 'iterations', 100, 'useGMRES', true);
end
figure;
tmp = cellfun(@(x) x.resvec, reports, 'uniformoutput', false);
tmp = [tmp{:}];
if size(tmp, 1) == 1
    tmp = [tmp; tmp];
end
names = cellfun(@(x) x.type, bases, 'uniformoutput', false);
semilogy(tmp, 'o-','MarkerFaceColor',[.8 .8 .8])
legend(names)
title('Convergence of GMRES')

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
