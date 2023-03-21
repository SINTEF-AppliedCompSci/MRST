function activeCells = cart2active(G, cartCells)
%Compute active cell numbers from linear Cartesian index.
%
% SYNOPSIS:
%   activeCells = cart2active(G, c)
%
% PARAMETERS:
%   G - Grid data structure as described by `grid_structure`.
%   c - List of linear Cartesian cell indices.
%
% RETURNS:
%   activeCells - Active cell numbers corresponding to the individual
%                 Cartesian cell numbers in `c`.
%
% NOTE:
%   This function provides the inverse mapping of the `G.cells.indexMap`
%   field in the grid data structure.
%
% SEE ALSO:
%   `grid_structure`, `gridLogicalIndices`

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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


assert (max(cartCells) <= prod(G.cartDims));
assert (min(cartCells) > 0);

actnum                   = false([prod(G.cartDims), 1]);
actnum(G.cells.indexMap) = true;

c           = cumsum(actnum) .* actnum;
activeCells = c(cartCells);
activeCells = activeCells(activeCells > 0);

%map = @(i) sparse(1:sum(a(i)), i(a(i)),1, sum(a(i)), numel(a))*cumsum(a);
end
