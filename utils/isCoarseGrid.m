function isCoarse = isCoarseGrid(G)
%Check if a grid is a coarse grid or a fine grid
%
% SYNOPSIS:
%   isCoarse = isCoarsegrid(G);
%
% REQUIRED PARAMETERS:
%   G        - Grid structure.
%
%
% RETURNS:
%   isCoarse - Boolean indicating of the grid is a coarse grid made from a
%              finer grid.
%
% SEE ALSO:
%   `plotGrid`, `generateCoarseGrid`

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
    isCoarse = isfield(G, 'parent') && isfield(G, 'partition');
end