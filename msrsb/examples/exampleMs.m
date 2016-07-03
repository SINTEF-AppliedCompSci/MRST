%% Simple conceptual multiscale example
% This example goes through the construction of a simple multiscale solver
% based on smoothed basis functions. We begin by defining the fine scale
% problem and plotting the permeability field.

% This value corresponds to the layers of the SPE10 dataset we are going to
% use. It can either be a single value or multiple values. Values from 1 to
% 35 will result in smooth, lognormal values while the values from 36 to 85
% will result in channelized formations that are considered to be
% significantly more challenging to upscale.

layerNo = 85;
% layerNo = 1;

mrstModule add mrst-gui spe10 coarsegrid incomp
[G, ~, rock] = SPE10_setup(layerNo);

% Boundary conditions corresponding to unit pressure drop over the entire
% domain
bc = [];
bc = pside(bc, G, 'ymin', 0);
bc = pside(bc, G, 'ymax', 1);

% Transmissibilities on fine scale
T = computeTrans(G, rock);
% Single phase unit fluid model. Since we are solving a single-phase
% incompressible problem without gravity, the values here does not matter.
fluid = initSingleFluid('rho', 1, 'mu', 1);
% Set up empty solution structure.
state = initResSol(G, 0);

% Plot the permeability field (log-scale).
figure(1); clf
plotToolbar(G, log10(rock.perm(:, 1)));
axis equal tight off
colormap jet
view(90, 90);

%% Solve fine scale problem
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