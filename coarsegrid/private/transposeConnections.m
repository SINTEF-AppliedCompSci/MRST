function [connPos, conns] = transposeConnections(N)
%Derive block->connection (face) mapping from connection->block mapping
%
% SYNOPSIS:
%   [connPos, conns] = transposeConnections(N)
%
% PARAMETERS:
%   N - Neighbourship definition.  An m-by-2 array representing connections
%       between entities.  May be a (coarse-scale) equivalent of the
%       faces.neighbors fields in the MRST grid structure.
%
% RETURNS:
%   connPos - 
%       Indirection array into the 'conns' array of block->connection
%       mappings.  Specifically, the connections of block 'i' are located
%       in elements
%
%              connPos(i) : connPos(i + 1) - 1
%
%       of the 'conns' array.
%
%   conns - 
%       Connection data ordered per block.
%
% NOTE:
%   If the connection/neighbourship definition 'N' indeed is the
%   faces.neighbors field of a grid_structure, then the (connPos,conns)
%   packed data array pair corresponds to (a possibly permuted version of)
%   the (cells.facePos,cells.faces) packed data array pair of the same grid
%   structure.
%
%   This function uses SORTROWS.
%
% SEE ALSO:
%   `coarseConnections`, `sortrows`, `grid_structure`.

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


   assert (isnumeric(N) && (size(N, 2) == 2), ...
          ['Neighbourship definition ''N'' must be an m-by-2 ', ...
           'numeric array of entity pairs.']);

   nconn = size(N, 1);

   V = sortrows([reshape(N, [], 1), repmat((1 : nconn) .', [2, 1])]);

   n = accumarray(V(:,1) + 1, 1);

   connPos = cumsum([1; n(2:end)]);
   conns   = V(n(1) + 1 : end, 2);
end
