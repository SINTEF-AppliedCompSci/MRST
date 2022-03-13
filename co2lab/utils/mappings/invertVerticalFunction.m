function h = invertVerticalFunction(f, g, z, y)
%Solve y(z)=f for piecewise linear, monotonically increasing function y.
%
% SYNOPSIS:
%   h = invertVerticalFunction(f, G, z, y)
%
% PARAMETERS:
%   f - Specific function values of the function y=y(z).  One scalar value
%       for each column in the top-surface grid.  As a special case, a
%       single, scalar value is repeated for each column in the grid.
%
%   G - A top-surface grid as defined by function 'topSurfaceGrid'.
%
%   z - Sample points (typically column depths) for the function y=y(z).
%       Possibly computed by means of function 'cumulativeHeight'.
%
%   y - Sampled values of the function y=y(z) at the points 'z'.  Often
%       corresponds to cumulative pore height as defined by function
%       'cumulativePoreHeight'.
%
% RETURNS:
%   h - Depth values (z) such that y(z)=f in each column.  One scalar value
%       for each column.
%
% EXAMPLE:
%   % Define geometry and rock data.
%   grdecl    = makeModel3([100, 60, 15]);
%   G         = computeGeometry(processGRDECL(grdecl));
%   G         = topSurfaceGrid(G);
%   rock.poro = rand([numel(G.columns.cells), 1]);        % Synthetic...
%
%   % Compute cumulative column depths and pore heights.
%   z  = cumulativeHeight(G);
%   ph = cumulativePoreHeight(G, rock);
%
%   % Compute depth (h) of 0.5 (total) column pore height.
%   f = 0.5 * ph(G.cells.columnPos(2 : end) - 1);
%   h = invertVerticalFunction(f, G, z, ph);
%
%   % Plot depth map of 0.5 pore height as a fraction of the
%   % corresponding (total) column depth.
%   plotCellData(G, h ./ z(G.cells.columnPos(2:end) - 1))
%   caxis([0, 1]); colorbar
%
% NOTE:
%   The function y=y(z) is assumed to be a piecewise linear, monotonically
%   increasing function (per column) satisfying y(0)=0.
%
% SEE ALSO:
%   `cumulativePoreHeight`, `cumulativeHeight`, `topSurfaceGrid`.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

% $Date: 2012-09-04 17:18:36 +0200 (Tue, 04 Sep 2012) $
% $Revision: 9535 $

   nc = reshape(diff(g.cells.columnPos), [], 1);
   if numel(f) == 1,
      % Scalar f => same value in all columns.
      f = f(ones(size(nc)));
   end

   % Exclude non-physical f<0 and guarantee col-vec shape.
   f = reshape(max(f, 0), [], 1);

   t      = rldecode([f, (1 : numel(f)).'], nc);
   t(:,1) = y - t(:,1);

   n       = accumarray(t(:,2), t(:,1) < 0, [numel(f), 1]);
   first   = n == 0 ;
   outside = n == nc;

   h = zeros(size(n));
   p = g.cells.columnPos(1 : end-1) + n;

   % 1) Solution in first cell.
   %    Function value y==f at fraction f/y of the interval [0, z(pf)).
   %
   pf         = p(first);
   h(first)   = (f(first) ./ y(pf)) .* z(pf); clear pf

   % 2) Solution outside (below) column.
   %    Arbitrarily assign maximum column depth.
   %    May want to issue a warning in this case...
   %
   h(outside) = z(p(outside) - 1);

   % 3) General case: Solution somewhere in the interior
   %    (neither first cell nor outside the column).
   %    Simply invert the linear function in this interval, [z(j-1), z(j)).
   %
   j    = ~(first | outside);  clear first outside
   pj   = p(j);
   h(j) = (z(pj-1).*t(pj,1) - z(pj).*t(pj-1,1)) ./ ...
          (         t(pj,1) -        t(pj-1,1));
end
