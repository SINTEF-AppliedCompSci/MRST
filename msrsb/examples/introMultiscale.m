%% Introduction to finite-volume multiscale methods
% We will consider a simple single-phase flow problem for a heterogeneous
% permeability distribution,
%
% $$ \nabla \cdot \vec{v} = q, \quad \vec{v} = - \mathbf{K} \nabla p$$
%
% In this example, we compare three different approaches for
% solving this system: 
% 
% * A fine-scale solution using a TPFA-type discretization
% * A solution based on a flow-based upscaling, which gives a solution on a
% coarse grid with the same type of TPFA-discretization.
% * A finite-volume multiscale solver, which gives both a coarse-scale and
% a fine-scale approximation, together with a fine-scale velocity field
% admissible for transport.
%
% This example uses several different modules.
mrstModule add coarsegrid spe10 msrsb incomp upscaling ad-core streamlines
%% Define grid and rock
% We take a small subset of the SPE 10, model 2 upscaling benchmark. We
% pick a small subset since we are interested in examining the fine-scale
% linear system, which quickly becomes large as more cells are included.
dims = [20, 20];
physDims = dims .*[20, 10]*ft;
% Create grid
G = cartGrid(dims, physDims);
G = computeGeometry(G);
% Take the rock structure from layer 15
rock = getSPE10rock(1:dims(1), 1:dims(2), 15);
rock.perm = rock.perm(:, 1);
hT = computeTrans(G, rock);
% Coarsen into blocks made up of 5 x 4 fine cells
cdims = [5, 4];
p = partitionUI(G, cdims);
% Alternatively, we could make an unstructured coarse grid using e.g. Metis
% p = partitionMETIS(G, hT, 20);

% Make coarse grid
CG = generateCoarseGrid(G, p);
% Add coarse geometry (centroids, normals, etc)
CG = coarsenGeometry(CG);

% Plot the permeability and coarse grid
figure('units', 'normalized', 'position', [0.1, .1, .8, 0.8]);
plotCellData(G, log10(rock.perm(:, 1)), 'edgec', 'none')
axis tight
plotGrid(CG,'FaceColor','none','LineWidth',2);
%% Define and solve the fine-scale flow problem
% Set up a pressure drop along the x-axis.
bc = pside([], G, 'West', 100*barsa);
bc = pside(bc, G, 'East',  50*barsa);
fluid = initSingleFluid('rho', 1, 'mu', 1);
state0 = initResSol(G, 0);
% Call solver, with matrix output in state enabled
state = incompTPFA(state0, G, hT, fluid, 'MatrixOutput', true, 'bc', bc);
n = G.cells.num;
%% Show the fine-scale discretization
% We plot the system matrix for the fine-scale system. We used a two-point
% flux approximation (TPFA) to solve the problem and consequently the
% matrix has a banded structure for this structured grid. As each cell in
% a Cartesian grid has four neighbors, we get four off-diagonals bands plus
% the main diagonal. As the incompressible pressure equation is linear, we
% can fully translate the governing equations into the form
% $$ A\mathbf{p} = \mathbf{q} $$
%
% In order to show this, we plot both the full system matrix and a
% zoomed-in view around the degrees of freedom corresponding to a cell in
% the middle of the domain.
drawSquare = @(x, y) patch(x([1, 1, 2, 2]), y([1, 2, 2, 1]), 0, ...
                            'EdgeColor', 'k', 'FaceColor', 'none');
% Get system matrix from state
A = full(abs(state.A));

dx = 1.5*dims(1);
dy = 1.5*dims(2);
% We inspect a cell in the middle of the domain
inspect = sub2ind(dims, ceil(dims(1)/2), ceil(dims(2)/2));
x = [inspect - dx, inspect + dx];
y = [inspect - dy, inspect + dy];
% Plot A
clf;
subplot(1, 2, 1)
imagesc(A);
colormap([1 1 1; parula(255)]);
axis([.5 n+.5 .5 n+.5]); axis square;
hold on;
drawSquare(x, y);
title('Fine-scale matrix')
% Plot zoomed view
subplot(1, 2, 2);
imagesc(A);
ylim(x);
xlim(y);
axis square
colormap([1 1 1; parula(255)]);
title('Local five-point structure');
%% Plot the fine-scale solution together with streamlines
% We plot the pressure field and trace five streamlines from each of the
% cells on the left boundary. Heterogeneities in the permeability field and
% variation in the underlying porosity results in deviation from linear
% flow, even with a linear pressure drop over the domain.
ijk = gridLogicalIndices(G);
cells = find(ijk{1} == 1);
p = repmat([0.25, 0.25; ...
            0.75, 0.25; ...
            0.25, 0.75; ...
            0.75, 0.75; ...
            0.5, 0.5], numel(cells), 1);
cno = rldecode(cells, 5);
startpos = [cno, p];

trace = @(state) pollock(G, state, startpos, 'pvol', rock.poro);

clf;
colormap default
plotCellData(G, state.pressure/barsa, 'edgec','none');
h = streamline(trace(state));
set(h, 'Color', 'w');
colormap(jet);
axis equal tight
title('Fine-scale solution');
colorbar('horiz')

%% Perform upscaling of the permeability field
% We create a upscaled permeability field with a flow-based upscaler from
% the 'upscale' module. In this routine, a series of local flow problem are
% used to get an effective coarse-scale permeability.
clf;
crock.perm = upscalePerm(G, CG, rock);

subplot(2, 1, 1);
Kf = log10(rock.perm(:,1));
plotCellData(G, Kf, 'edgec', 'none')
axis equal tight
c = caxis();
title('Fine permeability')

subplot(2, 1, 2);
Kc = log10(crock.perm(:,1));
plotCellData(CG, Kc, 'edgec', 'none')
axis equal tight
caxis(c);
title('Coarse permeability')
%% Solve the coarse problem
% The coarse problem is essentially the same governing equation as the
% fine-scale problem, but with an upscaled permeability field
% $\mathbf{K}_u$, and posed on the coarser mesh,
%
% $$ \nabla \cdot \vec{v}_u = q_u, \quad \vec{v}_u = - \mathbf{K}_u \nabla p_u$$
%
% With a linear system
% $$ A_u \mathbf{p}_u = \mathbf{q}_u $$
%
chT = computeTrans(CG, crock);
cbc = coarsenBC(CG, bc);
cstate = initResSol(CG, 0);
cstate = incompTPFA(cstate, CG, chT, fluid, 'MatrixOutput', true, 'bc', cbc);
%% Solve with a multiscale solver
% We now solve the problem with a multiscale solver. We do this by
% generating a set of basis functions $B$ from the system matrix. We can
% then solve the effective coarse scale problem to get the coarse
% multiscale pressure solution,
%
% $$ A_c \mathbf{p}_c = \mathbf{q}_c, \quad A_c = RAP, \, \mathbf{q}_c = R q $$
%
% where $R$ is the restriction operator and $B$ is a matrix representation
% of the basis functions. For a problem with $n$ fine cells and $n\times n$
% fine-scale system and a coarse grid with $m$ blocks, the basis functions
% and restriction operators take the form of rectangular matrices, 
%
% $$ B \in \mathbf{R}^n \times \mathbf{R}^m, \mathbf{R}^m \times {R}^n $$
%
% Due to the interpretation of the basis functions as interpolators from
% unit values, each entry of $B$ takes on values in $[0, 1]$. The
% restriction operator $R$ represents a integration of boundary fluxes and
% therefore the entries are either one or zero.
%
% We next generate basis functions, create a coarse system and visualize
% it. The intentionally small problem allows us to easily visualize the
% systems used to construct the coarse approximation. Note that we normally
% use functions such as 'incompMultiscale' to form this coarse system, but
% we explicitly write out the expressions here for pedagogical reasons.

basis = getMultiscaleBasis(CG, A);
Am = basis.R*state.A*basis.B;
bm = basis.R*state.rhs;
p_ms = Am\bm;
m = CG.cells.num;
clf;
subplot(4, 4, [5, 9, 13])

imagesc(full(basis.B));
colormap([1 1 1; parula(128)]);
daspect([1, 6, 1])
axis tight
title('P');

subplot(4, 4, 2:4);
imagesc(full(basis.R));
daspect([4, 1, 1]);
axis tight
title('R')

subplot(4, 4, [6:8, 10:12, 14:16]);
imagesc(full(abs(Am))); colormap([1 1 1; parula(255)]);
axis equal tight
th = text(dims(1), 3, 'A_c = RAP', 'FontSize', 18, 'HorizontalAlignment', 'right');
%% Compare the upscaled and multiscale coarse systems
% Fundamentally these two approaches are similar in that they, from
% fine-scale information, produce coarse scale discretizations which seek
% to account for the known fine-scale structure. We observe that the
% upscaled system has five major diagonals, just like the original system,
% while the multiscale solver has a more irregular pattern similar to
% what we would get if we discretized the coarse system with a multi-point
% scheme. As the system is formed numerically, all nine diagonals are not
% fully filled in, since the magnitude of the entries are dependent on the
% basis functions.
clf;
subplot(1, 2, 1);
imagesc(full(abs(cstate.A))); colormap([1 1 1; parula(255)]);
axis equal tight
title('Upscaled system');

subplot(1, 2, 2);
imagesc(full(abs(Am))); colormap([1 1 1; parula(255)]);
axis equal tight
title('Multiscale system');
%% Plot the solutions
% A key feature of multiscale methods is that we can create a fine-scale
% approximate pressure from the coarse system by applying the prolongation
% operator,
%  $$ \mathbf{p} \approx B \mathbf{p}_c $$
%
% We plot the fine-scale reference, the two coarse-scale approximations and
% the fine-scale approximation from the multiscale solver. Note that we use
% the same color axis for all pressures, as the coarse solvers generally do
% not have the outliers present in the fine-scale solution.
p_upscaled = cstate.pressure;
p_fine = state.pressure;
p_prolongated = basis.B*p_ms;
ca = [min(p_fine), max(p_fine)];

clf;
colormap default
subplot(2, 2, 1);
plotCellData(CG, p_upscaled);
axis equal tight;
caxis(ca);
title('Upscaled')

subplot(2, 2, 2);
plotCellData(G, p_fine);
axis equal tight;
caxis(ca);
title('Fine-scale')

subplot(2, 2, 3);
plotCellData(CG, p_ms);
axis equal tight;
caxis(ca);
title('Coarse MS');

subplot(2, 2, 4);
plotCellData(G, p_prolongated);
axis equal tight;
caxis(ca);
title('Fine MS');
%% Plot the error in the fine-scale multiscale solution
% The multiscale solution has two sources of error: The approximations made
% to the basis functions to ensure local support, and the error in the
% coarse scale system. In this case, however, the error is fairly low. We
% plot the scaled point-wise error
%  $$ e = \frac{|\mathbf{p} - \mathbf{p}_{ref}|}{\max p_{ref}} $$
clf;
plotCellData(G, abs(p_prolongated - p_fine)/max(p_fine));
axis equal tight;
caxis auto;
colorbar
title('Error, MS');
%% Reconstructed velocity field
% A fine-scale pressure approximation is generally of limited use for
% multiphase flow. We typically are interested in the solution of a
% transport equation,
% 
% $$ \frac{S \phi}{ \Delta t} + \Nabla \cdot (f(S) \vec{v}) = q
% 
% Under the limiting assumptions of incompressible, immiscible flow, we
% require a divergence free velocity field outside of source cells to solve
% this equation. The 'incompMultiscale' solver performs a flux
% reconstruction in order to produce a velocity field suitable for
% transport.
%
% We use the Divergence operator from 'ad-core' to verify this. If we plot
% the fine-scale divergence, we observe that the values are only nonzero
% near the boundary, where source terms are present. The same is the case
% for the multiscale solver when reconstruction is used, but it is not the
% case when we disable the flux reconstruction.
op = setupOperatorsTPFA(G, rock);

Div = @(flux) op.Div(flux(op.internalConn));
state_ms = incompMultiscale(state0, CG, hT, fluid, basis, 'bc', bc, 'reconstruct', true);
state_ms2 = incompMultiscale(state0, CG, hT, fluid, basis, 'bc', bc, 'reconstruct', false);

df = Div(state.flux);
dms = Div(state_ms.flux);
dms2 = Div(state_ms2.flux);

mv = max(max(abs(df)), max(abs(dms)));
ca = [-mv, mv];

cmap = interp1([0; 0.5; 1], [1, 0, 0; 1, 1, 1; 0, 0, 1], 0:0.01:1);

clf;
subplot(1, 3, 1)
plotCellData(G, dms);
colormap(cmap);
caxis(ca);
axis equal tight
colorbar('horiz')
title('Divergence (fine-scale)');

subplot(1, 3, 2);
plotCellData(G, dms);
colormap(cmap);
caxis(ca);
axis equal tight
colorbar('horiz')
title('Divergence (reconstructed velocity)')

subplot(1, 3, 3);
plotCellData(G, dms2);
colormap(cmap);
caxis(ca);
axis equal tight
colorbar('horiz')
title('Divergence (only multiscale)')
%% Plot fine-scale pressures side by side - together with streamlines
% We have excellent agreement between the fine-scale pressures and the
% velocity field when visualized as streamlines.
clf
colormap default
subplot(1, 2, 1);
plotCellData(G, p_fine);
axis equal tight;
h = streamline(trace(state));
set(h, 'Color', 'w');
title('Fine-scale');

subplot(1, 2, 2);
plotCellData(G, p_prolongated);
axis equal tight;
h = streamline(trace(state_ms));
set(h, 'Color', 'w');
title('Multiscale');