%% Introduction to multiscale finite-volume methods
% This example introduces you to the basic concepts of multiscale
% finite-volume methods. To this end, we consider a single-phase
% flow problem written in terms of conservation of mass and Darcy's law
%
% $$ \nabla \cdot \vec{v} = q, \quad \vec{v} = - \mathbf{K} \nabla p$$
%
% Alternatively, we could write this as a variable-coefficient Poisson
% equation,\nabla(K\nabla p)=q, which is one of the primary examples of a
% second-order elliptic partial differential equation. The multiscale
% behavior of this problem comes from the coefficient K, which represents
% permeability and usually has large local changes and spatial correlation
% structures that span a wide variety of scales.
%
% To discretize the problem, we use a standard finite-volume method with a
% two-point flux approximation. In a standard upscaling method, one first
% constructs a coarse grid and then computes effective permeability values
% that attempt to represent the flow equation in an averaged sense on this
% grid. With this, one can solve the flow problem approximately using fewer
% unknowns.
%
% In a multiscale method, one also tries to solve the flow problem using
% fewer unknowns, but instead of using effective coarse-scale
% permeabilities, we define a set of multiscale basis functions that
% represent how the pressure will distribute itself locally, given a unit
% pressure drop between the center of one coarse grid block and the centers
% of all the surrounding grid blocks. Using these basis functions, we can
% reduce the fine-scale discretization matrix to a much smaller matrix and
% solve this to construct an approximate solution on the coarse grid as
% well as on the fine grid. We can also reconstruct mass-conservative
% fluxes on the fine grid or on any partition thereof.
%
% The example compares three different approximate solutions to the flow
% problem just described:
% 
% * A fine-scale solution 
% * A coarse-scale solution computed using flow-based upscaling
% * A multiscale solution computed with the MsRSB multiscale method

% This example uses several different modules.
mrstModule add coarsegrid spe10 msrsb incomp upscaling ad-core streamlines

%% Define grid and rock
% As an example of a strongly heterogeneous problem, we consider
% permeability sampled from Model 2 of the SPE10 upscaling benchmark. We
% pick a small subset, since we are interested in examining the fine-scale
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
% Alternatively, we could make an unstructured coarse grid using, e.g., Metis
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
title('Permeability (log-scale) and coarse partition');

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
% flux approximation (TPFA) to solve the problem, and consequently the
% matrix has a banded structure for this structured grid. As each cell in
% a Cartesian grid has four neighbors, we get four off-diagonals bands plus
% the main diagonal. The incompressible pressure equation is linear, and we
% can thus fully translate the governing equations into the form
% $$ A\mathbf{p} = \mathbf{q} $$
%
% To show this, we plot both the full system matrix and a zoomed-in view
% around the degrees of freedom corresponding to a cell in the middle of
% the domain.
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
% variation in the underlying porosity result in deviation from linear
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
% We create an upscaled permeability field with a flow-based upscaler from
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
% The coarse problem consists of the same governing equation as the
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
% We now solve the problem with a multiscale solver. We do this by first
% generating a set of basis functions $B$ from the system matrix. We can
% then form and solve a reduced problem to get the coarse multiscale
% pressure solution,
%
% $$ A_c \mathbf{p}_c = \mathbf{q}_c, \quad A_c = RAP, \, \mathbf{q}_c = R q $$
%
% where $R$ is the restriction operator and $B$ is a matrix representation
% of the basis functions. For a problem with $n$ fine cells, giving an
% $n\times n$ fine-scale system, and a coarse grid with $m$ blocks, the
% basis functions and restriction operators take the form of rectangular
% matrices,
%
% $$ B \in \mathbf{R}^n \times \mathbf{R}^m, \mathbf{R}^m \times {R}^n $$
%
% Due to the interpretation of the basis functions as interpolators from
% unit values, each entry of $B$ takes values in $[0, 1]$. The
% restriction operator can be defined in two different ways. First of all,
% we could set R as the transpose of P, which would give us a Galerkin
% formulation. Herein, we will instead use a finite-volume restriction that
% simply sums the rows corresponding to all cells that make up each coarse
% block. Thus R takes values in {0,1}.
%
% We next generate basis functions, create a coarse system, and visualize
% it. Since the problem is very small, we can easily visualize the systems
% used to construct the coarse approximation. Note that we normally use
% functions such as 'incompMultiscale' to form this coarse system, but here
% we explicitly write out the expressions for pedagogical reasons.

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
% Fundamentally, these two approaches are similar in that they, from
% fine-scale information, produce coarse-scale discretizations that seek to
% account for the known fine-scale structure. We observe that the upscaled
% system has five major diagonals, just like the original system, whereas
% the multiscale solver has a more irregular pattern similar to what we
% would get if we discretized the coarse system with a multipoint scheme.
% As the system is formed numerically, all nine diagonals are not fully
% filled in, since the magnitude of the entries depend on the basis
% functions.
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
% We plot the fine-scale reference, the two coarse-scale approximations,
% and the fine-scale approximation from the multiscale solver. We enforce
% the same color axis for all subplots, as the coarse solvers generally do
% not have the outliers present in the fine-scale solution.
p_upscaled = cstate.pressure;
p_fine = state.pressure;
p_prolongated = basis.B*p_ms;
ca = [min(p_fine), max(p_fine)]/barsa;

clf;
colormap default
subplot(2, 2, 1);
plotCellData(CG, p_upscaled/barsa,'EdgeColor','none'); 
plotFaces(CG,1:CG.faces.num)
axis equal tight;
caxis(ca);
title('Upscaled pressure')

subplot(2, 2, 2);
plotCellData(G, p_fine/barsa);
axis equal tight;
caxis(ca);
title('Fine-scale pressure')

subplot(2, 2, 3);
plotCellData(CG, p_ms/barsa,'EdgeColor','none'); 
plotFaces(CG,1:CG.faces.num)
axis equal tight;
caxis(ca);
title('Coarse MS pressure');

subplot(2, 2, 4);
plotCellData(G, p_prolongated/barsa);
axis equal tight;
caxis(ca);
title('Fine MS pressure');

%% Plot the error in the fine-scale multiscale solution
% The multiscale solution has two sources of error: The approximations made
% to ensure that the basis functions have local support, and the error in
% the coarse-scale system. In this case, however, the error is fairly low.
% We plot the scaled point-wise error
%  $$ e = \frac{|\mathbf{p} - \mathbf{p}_{ref}|}{\max p_{ref}} $$
clf;
plotCellData(G, abs(p_prolongated - p_fine)/max(p_fine));
axis equal tight;
caxis auto;
colorbar
title('Relative MS pressure error');

%% Reconstructed velocity field
% A fine-scale pressure approximation is usually not our primary interest
% when studying a multiphase flow problem. Instead, we are typically
% interested in the solution of a transport equation,
% 
% $$ \frac{S \phi}{ \Delta t} + \Nabla \cdot (f(S) \vec{v}) = q
% 
% Under the limiting assumptions of incompressible, immiscible flow, we
% require a divergence free velocity field outside of source cells to solve
% this equation. The 'incompMultiscale' solver performs a flux
% reconstruction in order to produce a velocity field suitable for
% transport.
%
% We use the discrete divergence operator from 'ad-core' to verify this. If
% we plot the fine-scale divergence, we observe that the values are only
% nonzero near the boundary, where source terms are present. The same is
% the case for the multiscale solver when reconstruction is used, but it is
% not the case when we disable the flux reconstruction.
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
plotCellData(G, p_fine,'EdgeColor','none');
axis equal tight;
h = streamline(trace(state));
set(h, 'Color', 'k','linewidth',.5);
title('Fine-scale');

subplot(1, 2, 2);
plotCellData(G, p_prolongated,'EdgeColor','none');
axis equal tight;
h = streamline(trace(state_ms));
set(h, 'Color', 'k','linewidth',.5);
title('Multiscale');

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
