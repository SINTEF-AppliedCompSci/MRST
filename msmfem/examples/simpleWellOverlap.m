%% Multiscale Pressure Solver with Overlap:
% Use overlap in multiscale basis function to improve accuracy of the
% simulation.  We compare the fine-grid pressure solver and multiscale
% pressure solvers with and without overlap by solving the single-phase
% pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\frac{K}{\mu}\nabla p,$$
%
% for a Cartesian grid with anisotropic, homogeneous permeability.
%
% This example shows the use of overlap around the well.

mrstModule add coarsegrid mimetic msmfem incomp
%% What is overlap?
% Overlap means that we extend the support of the multiscale basis
% functions, that is, we include more fine cells when computing the basis
% function. Using overlapping domains enables capturing more complex flow
% patterns, particularly for very coarse grids, at the expense of increased
% coupling in the resulting systems of linear equations. The use of overlap
% is particularly attractive in cases of highly heterogeneous reservoirs of
% cases where wells are placed near the corner of a coarse well-block. In
% the following, we describe the different types of overlap implemented in
% the Matlab Reservoir Simulation Toolbox. Let the thick lines be the
% coarse grid, the thin lines represent the fine grid, the shaded region
% the overlap (of fine cells), and the red circle be the well.
%
% Overlap in the reservoir basis function for a coarse face ij means that
% we include fine cells around the the neighboring coarse blocks when we
% compute the basis for face ij:
%
% <<basisIJ.png>>
%
% The flow near a well is described by a multiscale well basis. This
% basis normally only has support in the coarse well-block. If the well
% lies close to the corner of the well-block, this will be a poor
% description of the flow near the well because it will not describe the
% flow over the corner. In such cases we can use overlap around
% the well ('overlapWell'), here shown in different colors with increasing
% number of overlap:
%
% <<basisOverlapWell.png>>
%
% In some cases (e.g in case of highly hetereogeneous permeability) it
% might be necessary to use overlap around the entire well block
% ('OverlapBlock'):
%
% <<basisOverlapBlock.png>>
%
%

%% Define and visualize the model
% We set up a simple 40-by-40-by-2 cartesian grid model with anisotropic,
% homogeneous permeability and assume that the reservoir is filled with a
% single fluid (with default values).
G         = cartGrid([40, 40, 2]);
G         = computeGeometry(G);
diag_k    = [500, 500, 50] .* milli*darcy;
rock.perm = repmat(diag_k, [G.cells.num, 1]);
fluid     = initSingleFluid('mu' ,    1*centi*poise     , ...
                            'rho', 1014*kilogram/meter^3);

gravity off

% Set two vertical wells.
W = verticalWell([], G, rock,  8,  8, 1:G.cartDims(3), ...
                 'InnerProduct', 'ip_simple', ...
                 'Type', 'rate', 'Val', 1*meter^3/day, ...
                 'Radius', 0.1, 'Comp_i', 1);
W = verticalWell(W , G, rock, 33, 33, 1:G.cartDims(3), ...
                 'InnerProduct', 'ip_simple', ...
                 'Type', 'bhp', 'Val', 0, 'Radius', 0.1, 'Comp_i', 0);

%% Partition the grid and compute mimetic innerproduct
% We partition the fine grid uniformly into 5-by-5-by-1 blocks using no
% overlap in the basis reservoir basis functions.
part = partitionCartGrid(G.cartDims, [5, 5, 1]);
CG   = generateCoarseGrid(G, part);
S    = computeMimeticIP(G, rock);

% Note: Only use overlap around wells, not in reservoir basisfunctions
mu  = fluid.properties(initResSol(G, 0.0));
kr  = fluid.relperm(ones([G.cells.num, 1]), initResSol(G, 0.0));
mob = kr ./ mu;
CS  = generateCoarseSystem(G, rock, S, CG, mob, 'Overlap', 0);

clf,
   plotCellData(G, mod(part,2), 'EdgeColor', 'k');
   outlineCoarseGrid(G, part, 'LineWidth', 3);
   plotWell(G, W, 'radius', 0.1, 'height', 2);
   view(3), axis tight off

%% Generate coarse well system with overlap around the well
% As can be seen from the plot of the model, the wells are near a corner of
% their supporting coarse blocks, which may lead to highly inaccurate
% results (as we will see shortly). To impreove the accuracy, we can
% introduce onverlap in the well basis functions. As pointed out above, we
% have two different choices: 'OverlapWell' and 'OverlapBlock'. The first
% includes cells around the well and the latter includes cells around the
% well-block (the coarse block containing the well). We apply overlap
% around the well for 'W1' and use no overlap for 'W' as a reference.
%
W1 = W;

% First no overlap
W  = generateCoarseWellSystem(G, S, CG, CS, mob, rock, W, ...
                              'OverlapWell',  0, 'OverlapBlock', 0);
% Then overlap around well
W1 = generateCoarseWellSystem(G, S, CG, CS, mob, rock, W1, ...
                              'OverlapWell', 10, 'OverlapBlock', 0);

%% Solve the global flow problems
%  The mass matrix B in the linear system will not be block-diagonal when
%  we use overlap. Thus, the lineary system system can not be reduced by
%  Schur-complement reduction as we do in the hybrid solver. Therefore, the
%  coarse system with overlap must be solved with option 'Solver', 'mixed'
%  as input to solveIncompFlowMS.

% Fine scale reference solution:
ref = initState(G, W, 0, 1);
ref = incompMimetic(ref, G, S, fluid, 'wells', W);

% Coarse scale - no overlap:
s = initState(G, W, 0, 1);
s = solveIncompFlowMS(s,  G, CG, part, S, CS, fluid, ...
                      'wells', W, 'Solver', 'mixed');

% Coarse scale - overlap:
s1 = initState(G, W, 0, 1);
s1 = solveIncompFlowMS(s1, G, CG, part, S, CS, fluid, ...
                       'wells', W1, 'Solver', 'mixed');

% Report pressure drop between wells.
dp = @(x) num2str(convertTo(x.wellSol(1).pressure, barsa()));
disp(['DeltaP - Fine:          ', dp(ref)])
disp(['DeltaP - Ms no overlap: ', dp(s  )])
disp(['DeltaP - Ms overlap:    ', dp(s1 )])


%% Plot solution
cellNo    = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';
diverg    = @(v) accumarray(cellNo, abs(convertTo(v, meter^3/day)));
plot_flux = @(x) plotCellData(G, diverg(x), 'EdgeAlpha', 0.125);
well_bas  = @(w) S.BI * sparse(w.CS.basis{1}{1}, 1, w.CS.basis{1}{3}, ...
                               size(G.cells.faces,1), 1);

subplot(2,3,1)
   plot_flux(ref.flux(G.cells.faces(:,1)));
   title('Flux - fine scale'); cx = caxis;

subplot(2,3,2)
   plot_flux(s.flux(G.cells.faces(:,1)))  ;
   title('Flux - no overlap'); caxis(cx)

subplot(2,3,3)
   plot_flux(s1.flux(G.cells.faces(:,1))) ;
   title('Flux - overlap')   ; caxis(cx)

subplot(2,3,5)
   plot_flux(full(well_bas(W(1))));
   axis([0, 20, 0, 20]);  cx = caxis;
   title('Well basis - no overlap')

subplot(2,3,6)
   plot_flux(full(well_bas(W1(1))));
   axis([0, 20, 0, 20]); caxis(cx);
   title('Well basis - overlap')

%%
displayEndOfDemoMessage(mfilename)

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
