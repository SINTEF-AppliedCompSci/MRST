%% Simple Conceptual Multiscale Example (Layers of SPE10)
% This example goes through the construction of a simple multiscale solver
% based on restriction-smoothed basis functions. We begin by defining the
% fine-scale problem and plotting the permeability field, which represents
% a number of layers sampled from the SPE10 model. Layers from 1 to 35 will
% result in smooth, lognormal permeabilities, whereas layers 36 to 85 will
% result in channelized formations that are considered to be significantly
% more challenging for upscaling and multiscale methods.
%
% We solve the flow problem using three methods: a fine-scale TPFA solver
% from the 'incomp' module, the multiscale finite-volume (MsFV) method from
% the 'msfv' module, and the MsRSB method implemented in this module.

%%  Set up the fine-scale model
layerNo = 85;
% layerNo = 1;

mrstModule add mrst-gui spe10 coarsegrid incomp
[G, ~, rock] = getSPE10setup(layerNo);

% Transmissibilities on fine scale
T = computeTrans(G, rock);

% Boundary conditions corresponding to unit pressure drop over the entire
% domain
bc = [];
bc = pside(bc, G, 'ymin', 0);
bc = pside(bc, G, 'ymax', 1);

% Single-phase unit fluid model. Since we are solving a single-phase
% incompressible problem without gravity, the values here do not matter.
fluid = initSingleFluid('rho', 1, 'mu', 1);

% Set up empty solution structure.
state = initResSol(G, 0);

% Plot the permeability field (log-scale).
figure(1); clf
plotToolbar(G, log10(rock.perm(:, 1)));
axis equal tight off
colormap jet
view(90, 90);

%% Solve fine-scale problem
state = incompTPFA(state, G, T, fluid, 'MatrixOutput', true, 'bc', bc);

%% Set up coarsegrid
mrstModule add msrsb

coarsen = [5 10 5];
coarsedims = ceil(G.cartDims./coarsen);

% Generate partition vector
p = partitionUI(G, coarsedims);
% Grid structure from partition vector
CG = generateCoarseGrid(G, p);
% Add centroids / geometry information on coarse grid
CG = coarsenGeometry(CG);
% Store the support of each cell (used by multiscale basis construction)
CG = storeInteractionRegionCart(CG);

mrstModule add msfvm
DG = partitionUIdual(CG, coarsedims);
DG = makeExplicitDual(CG, DG);
CG.dual = DG;

%% Set up basis functions
A = getIncomp1PhMatrix(G, T);
% B = iteratedJacobiBasis(A, CG);
% R = controlVolumeRestriction(CG.partition);
% basis = struct('B', B, 'R', R);
basis_sb = getMultiscaleBasis(CG, A, 'type', 'rsb');
basis_fv = getMultiscaleBasis(CG, A, 'type', 'msfv');

basises = {basis_sb, basis_fv};

%% Set up smoother function
fn = getSmootherFunction('type', 'ilu');

%% Compute multiscale solutions
nb = numel(basises);

[states, reports] = deal(cell(nb, 1));
    
for i = 1:nb
    basis = basises{i};
    states{i} = incompMultiscale(state, CG, T, fluid, basis, 'bc', bc);
end

%% Solve using MS-GMRES
for i = 1:nb
    basis = basises{i};
    [~, reports{i}] = incompMultiscale(state, CG, T, fluid, basis, 'bc', bc,...
        'getSmoother', fn, 'iterations', 100, 'useGMRES', true);
end


%% Plot results
close all
figure;
plotToolbar(G, state.pressure)
colormap jet
view(90, 90)
axis tight off
title('Fine scale')
for i = 1:nb
    figure;
    plotToolbar(G, states{i}.pressure)
    colormap jet
    view(90, 90)
    axis tight off
    caxis([min(state.pressure), max(state.pressure)])
    title(basises{i}.type)
end

figure;
tmp = cellfun(@(x) x.resvec, reports, 'uniformoutput', false);
tmp = [tmp{:}];
if size(tmp, 1) == 1
    tmp = [tmp; tmp];
end
names = cellfun(@(x) x.type, basises, 'uniformoutput', false);
semilogy(tmp, '*-')
legend(names)
title('Convergence of GMRES')

%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2018 SINTEF ICT, Applied Mathematics.
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
