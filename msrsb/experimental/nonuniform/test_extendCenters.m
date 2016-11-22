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
[G, ~, rock] = getSPE10setup(layerNo);
G = cartGrid(G.cartDims(1:2));
G = computeGeometry(G);
rock.perm = rock.perm(:, 1);

rock = makeRock(G, 1, 1);

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
ref = incompTPFA(state, G, T, fluid, 'MatrixOutput', true, 'bc', bc);

%% Set up coarsegrid
mrstModule add new-multiscale

% coarsen = [5 10 5];
coarsen = [10, 20, 1];
coarsedims = ceil(G.cartDims./coarsen(1:G.griddim));

% Generate partition vector
% p = partitionUI(G, coarsedims);

p = partitionUniformPadded(G, coarsedims);
% Grid structure from partition vector
CG = generateCoarseGrid(G, p);
% Add centroids / geometry information on coarse grid
CG = coarsenGeometry(CG);
% Store the support of each cell (used by multiscale basis construction)
CG = addCoarseCenterPoints(CG);

% CG = storeInteractionRegion(CG, 'localTriangulation', true);

CG = storeInteractionRegionCart(CG);
%% Set up basis functions
A = getIncomp1PhMatrix(G, T);
basis = getMultiscaleBasis(CG, A, 'type', 'rsb');

ms = incompMultiscale(state, CG, T, fluid, basis, 'bc', bc);


%%
norm(ms.pressure - ref.pressure, 2)/norm(ref.pressure, 2)
norm(ms.pressure - ref.pressure, inf)/norm(ref.pressure, inf)

%% Plot results
close all
figure;
plotToolbar(G, ref.pressure)
colormap jet
view(90, 90)
axis tight off
title('Fine scale')

figure;
plotToolbar(G, ms.pressure)
colormap jet
view(90, 90)
axis tight off
caxis([min(ref.pressure), max(ref.pressure)])
title('Multiscale')
outlineCoarseGrid(G, CG.partition)
%%
close all
plotToolbar(G, basis);
axis tight off
view(-90, 90)
outlineCoarseGrid(G, CG.partition);
plotGrid(G, CG.cells.centers)
h = nan
for i = 1:CG.cells.num
    if ishandle(h)
        delete(h);
    end
    h = plotGrid(G, CG.cells.interaction{i});
    pause()
end
