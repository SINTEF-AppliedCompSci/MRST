function plotContours(g, value, n)
%Plot contours of cell data.
%
% SYNOPSIS:
%   plotContours(G, field, n)
%
% DESCRIPTION:
%   Function 'plotContours' is a poor-man's implementation of the
%   contour-level drawing function `contourf` designed to work with MRST's
%   grid structure and cell-based values rather than the pure Cartesian
%   (tensor product) nodal values of `contourf`.  The contour lines follow
%   the grid lines and are, typically, not smooth.
%
%   The MRST function `plotGridVolumes`, which interpolates the data field
%   onto a structured grid, is an alternative option for visualising levels
%   and iso-surfaces.
%
% PARAMETERS:
%   G     - Grid structure.
%
%   field - Scalar field (e.g., cell pressures).  Used as value data for
%           determining the contour locations.  One scalar value for each
%           cell in 'G'.
%
%   n     - Number of (equidistant) contour lines/levels for `field`.
%           Positive integer.
%
% RETURNS:
%   Nothing.
%
% EXAMPLE:
%   % Visualise the iso-levels of the pressure field of a quarter five-spot
%   % configuration.  Note that the actual pressure values in this case are
%   % artificial due to (very) high permeabilities (approximately 1e+12 D).
%   %
%   G     = computeGeometry(cartGrid([50, 50, 1]));
%   rock  = struct('perm', ones([G.cells.num, 1]));
%   T     = computeTrans(G, rock);
%   fluid = initSingleFluid('mu', 1, 'rho', 0);
%   src   = addSource([], [1, G.cells.num], [1, -1]);
%   x     = initState(G, [], 0);
%   x     = incompTPFA(x, G, T, fluid, 'src', src);
%
%   plotContours(G, x.pressure, 20), axis equal tight
%
% NOTE:
%   Function `plotContours` is only supported in three space dimensions
%   i.e., if `G.griddim == 3`.
%
% SEE ALSO:
%   `plotGridVolumes`, `plotCellData`, `plotFaces`, `contourf`

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


   assert (g.griddim == 3, ...
          ['Function ''%s'' is only supported in ', ...
           'three space dimensions.'], mfilename);

   levels = linspace(min(value), max(value), n + 1);
   cmap   = colormap;
   cmapix = fix(linspace(1, size(cmap,1), n));

   plotf = @(f, c) ...
      plotFaces(g, f, c, 'EdgeColor', 'none', ...
                'FaceAlpha', 0.4, 'outline', true);

   ext          = any(g.faces.neighbors == 0, 2);
   contourfaces = @(I) any(I(g.faces.neighbors + 1), 2);

   for i = 1 : n
      % Fill space between contours.
      I = [false; value >= levels(i) & value <= levels(i + 1)];
      f = find(contourfaces(I) & ext);

      colour = cmap(cmapix(i), :);

      plotf(f, colour); %#ok

      % Plot internal level sets.
      I = [false; value < levels(i + 1)];
      J = [false; value > levels(i + 1)];
      f = find(contourfaces(I) & contourfaces(J));

      plotf(f, colour); %#ok
   end
end
