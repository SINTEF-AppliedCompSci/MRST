function conn = blockConnectivity(N)
%Build block-to-neighbours map by transposing neighbourship definition
%
% SYNOPSIS:
%   conn = blockConnectivity(N)
%
% PARAMETERS:
%   N - Neighbourship definition.  An m-by-2 array of explicit connections,
%       akin to the 'faces.neighbors' field of the grid structure.
%
% RETURNS:
%   conn - Connectivity structure.  A B-by-1 cell array of
%          block-to-neighbour associations.  Specifically, the neighbours
%          of block 'b' are found in conn{b}.
%
% SEE ALSO:
%   `blockNeighbours`, `blockNeighbourship`.

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


   conn = accumarray(reshape(N        , [], 1),     ...
                     reshape(fliplr(N), [], 1), [], ...
                     @(x) { x });
end
