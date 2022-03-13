function [f,df] = integrateVertically(fun, h, g)
%Compute integral of function, vertically per column.
%
% SYNOPSIS:
%   F = integrateVertically(f, h, G)
%
% DESCRIPTION:
%   Computes a first order approximation to the integral
%
%        F(x,y) = \int_0^h f(x,y,z) dz
%
%   within each column of a top-surface grid.
%
% PARAMETERS:
%   fun - Sampled values of function f(x,y,z).  One value for each cell in
%         the underlying three-dimensional model.
%
%   h   - Upper bound for the integral (*).  Typically corresponds to the
%         surface depth.  One scalar value for each column in the
%         top-surface grid.  As a special case a single, scalar value is
%         repeated for each column in the top-surface grid.
%
%   G   - A top-surface grid as defined by function 'topSurfaceGrid'.
%
% RETURNS:
%   F - Vertically averaged function value, F(x,y).  One scalar value for
%       each column in the top-surface grid 'G'.
%
% SEE ALSO:
%   `topSurfaceGrid`, `cumulativeHeight`.

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

% $Date: 2012-01-30 11:39:51 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9019 $

   %fun = reshape(fun, size(g.columns.dz));
   [n, t] = fillDegree(h, g);

   % 1) Completely filled cells.
   p    = g.cells.columnPos(1 : end-1);
   i    = n > 0;
   j    = mcolon(p(i), p(i) + n(i) - 1);
   subs = rldecode(find(i), n(i));
   vals = fun(g.columns.cells(j)) .* g.columns.dz(j);

   % 2) Partially filled cells (one per column).
   %    Exclude completely filled columns which are already handled in the
   %    above case (case 1).
   %
   p    = p + n;
   i    = p < g.cells.columnPos(2 : end);
   subs = [subs; find(i)];
   vals = [vals; t(i) .* fun(g.columns.cells(p(i))) .* g.columns.dz(p(i))];

   % Evaluate integral in all columns.
   f = accumarray(subs, vals, size(n));

   % Derivative, either a partially filled cell or bottom cell in column
   ix = g.cells.columnPos(2 : end)-1;
   p(~i) = ix(~i);
   df = fun(g.columns.cells(p));
end
