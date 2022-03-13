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
% Overlap means that we extend the support of the multiscale basis
% functions, that is, we include more fine cells when computing the basis
% function. The use of overlap is particularly attractive in cases of
% highly heterogeneous reservoirs and wells that are placed near the corner
% of a coarse well-block. In the following, let the dashed lines be the
% coarse grid, the dotted lines represent the overlap, and let '*'
% represent a well.
%
%  Overlap in the reservoir basis function for a coarse face ij means that
% we include fine cells around the the neighboring coarse blocks when we
% compute the basis for face ij:
%
%                  . . . . . . . . . .
%                . o-----------------o .
%                . |        |        | .
%                . |   i    |   j    | .
%                . |        |        | .
%                . o-----------------o .
%                  . . . . . . . . . .
%
%  The flow near a well is described by a multiscale well basis. This
% basis normally only has support in the coarse well-block. If the well
% lies close to the corner of the well-block, this will be a poor
% description of the flow near the well because it will not describe the
% flow over the corner. In such cases we can use overlap around
% the well ('overlapWell'):
%
%                           . . . . .
%                           .       .
%                      o--------o   .
%                      |       *|   .
%                      |        | . .
%                      |        |
%                      o--------o
%
%  In some cases (e.g in case of highly hetereogeneous permeability) it
% might be necessary to use overlap around the entire well block
% ('OverlapBlock'):
%
%                    . . . . . . . .
%                    . o---------o .
%                    . |         | .
%                    . |    *    | .
%                    . |         | .
%                    . o---------o .
%                    . . . . . . . .
%
%  This example shows the use of overlap around the well.

mrstModule add coarsegrid mimetic msmfem incomp
%% Define and visualize the model
nx = 40; ny = 40; nz = 2;
cellDims = [nx, ny, nz];

G = cartGrid(cellDims);
G = computeGeometry(G);
diag_k     = [500, 500, 50] .* (1e-3*darcy());
rock.perm  = repmat(diag_k, [G.cells.num, 1]);
fluid = initSingleFluid('mu' ,    1*centi*poise     , ...
                        'rho', 1000*kilogram/meter^3);

% Set two vertical wells.
injCells  = [ 8,  8];
prodCells = [33, 33];

W = struct([]);
W = verticalWell(W, G, rock, injCells(1), injCells(2), 1:2, ...
                 'InnerProduct', 'ip_simple', ...
                 'Type', 'rate', 'Val', 1*meter^3/day,      ...
                 'Radius', 0.1);
W = verticalWell(W, G, rock, prodCells(1), prodCells(2), 1:2, ...
                 'InnerProduct', 'ip_simple', ...
                 'Type', 'bhp', 'Val', 0, 'Radius', 0.1);

% Visualize the model.
clf
plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [0.65, 0.65, 0.65]);
plotWell(G, W, 'radius', 0.1, 'color', 'r');
view(3), axis equal tight off



%% Partition the grid and assemble linear systems
part  = partitionCartGrid(cellDims, [5, 5, 1]);
CG    = generateCoarseGrid(G, part, 'Verbose', true);

gravity off
S  = computeMimeticIP(G, rock, 'Verbose', true);

% Note: Only use overlap around wells, not in reservoir basisfunctions
CS = generateCoarseSystem (G, rock, S, CG, ones([G.cells.num,1]), ...
                           'Verbose', true, 'Overlap', 0);

%% Generate coarse well system with overlap around the well
% For overlap in the well basis we has two different choices: 'OverlapWell'
% and 'OverlapBlock'. The first includes cells around the well and the
% latter cells around the well-block (the coarse block containing the
% well). We apply overlap around the well for 'W1' and use no overlap for
% 'W' as a reference.
W1 = W;

% Calculate total mobility (needed by generateCoarse*System).
state = initResSol(G, 0, 1);
mu    = fluid.properties(state);
s     = fluid.saturation(state);
kr    = fluid.relperm(s, state);

totmob = sum(bsxfun(@rdivide, kr, mu), 2);    clear mu s kr state

% First no overlap
W  = generateCoarseWellSystem(G, S, CG, CS, totmob, rock, W, ...
                              'OverlapWell',  0, 'OverlapBlock', 0);
% Then overlap around well
W1 = generateCoarseWellSystem(G, S, CG, CS, totmob, rock, W1, ...
                              'OverlapWell', 10, 'OverlapBlock', 0);

%% Solve the global flow problems
%  The mass matrix B in the linear system will not be block-diagonal when
%  we use overlap. Thus, the lineary system system can not be reduced by
%  Schur-complement reduction as we do in the hybrid solver. Therefore, the
%  coarse system with overlap must be solved using option 'Solver', 'mixed'
%  as input to solveIncompFlowMS.

% Fine scale reference solution:
xRef = initState(G, W, 0);
xRef = incompMimetic(xRef, G, S, fluid, 'wells', W);

% Coarse scale - no overlap:
xMs   = initState(G, W, 0);
xMs   = solveIncompFlowMS(xMs, G, CG, part, S, CS, fluid, ...
                          'wells', W, 'Solver', 'mixed');

% Coarse scale - overlap:
xMs1  = initState(G, W, 0);
xMs1  = solveIncompFlowMS(xMs1, G, CG, part, S, CS, fluid, ...
                          'wells', W1, 'Solver', 'mixed');

% Report pressure drop between wells.
dp = @(x) num2str(convertTo(x.wellSol(1).pressure, barsa()));
disp(['DeltaP - Fine:          ', dp(xRef)])
disp(['DeltaP - Ms no overlap: ', dp(xMs )])
disp(['DeltaP - Ms overlap:    ', dp(xMs1)])


%% Plot solution
cellNo    = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';
diverg    = @(v) accumarray(cellNo, abs(convertTo(v, meter^3/day)));
plot_flux = @(v) plotCellData(G, diverg(v), 'EdgeAlpha', 0.125);
well_bas  = @(w) S.BI * sparse(w.CS.basis{1}{1}, 1, w.CS.basis{1}{3}, ...
                               size(G.cells.faces,1), 1);
figure;
subplot(2,3,1)
   plot_flux(faceFlux2cellFlux(G, xRef.flux));
   title('Flux - fine scale'); cx = caxis;

subplot(2,3,2)
   plot_flux(faceFlux2cellFlux(G, xMs.flux ));
   title('Flux - no overlap'); caxis(cx)

subplot(2,3,3)
   plot_flux(faceFlux2cellFlux(G, xMs1.flux));
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
