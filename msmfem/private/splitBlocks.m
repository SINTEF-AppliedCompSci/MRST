function p = splitBlocks(blocks, p, splitFact, G)
%Split a set of coarse blocks into sub-blocks.
%
% SYNOPSIS:
%   p = splitBlocks(blocks, p, splitFact, G)
%
% PARAMETERS:
%   blocks    - Vector of coarse block numbers.  Each of these blocks is
%               split according to 'splitFact'.
%
%   p         - Original fine grid partition vector as defined by (e.g.)
%               partitionCartGrid or partitionUI.
%
%   splitFact - Vector of coarse block splitting factors in each physical
%               dimension.  Each coarse block in 'blocks' is split into
%               PROD(splitFact) new coarse blocks.
%
%   G         - Grid structure as described in grid_structure.
%
% RETURNS:
%   p - Updated partition vector such that all fine cells within the
%       original 'blocks' get (potentially) new block numbers.
%
% SEE ALSO:
%   `partitionCartGrid`, `refinePartitionForWells`, `compressPartition`.

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


if ~(isfield(G.cells, 'ijkMap') || ...
     (isfield(G, 'cartDims') && strcmp(G.type, 'tensorgrid'))),
   error(id('cellIJK:NotPresent'), ...
         'splitBlocks only supported for logically Cartesian grids.')
end

cIJK = getCellIJK(G);
nb   = max(p);

blocks    = reshape(blocks, 1, []);
splitFact = reshape(splitFact, 1, []);

for b = blocks,
   ix   = find(p == b);
   locp = splitBlock(cIJK(ix, :), splitFact) - 1;

   p(ix(locp > 0)) = nb + locp(locp > 0);

   nb = nb + max(locp);
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function p = splitBlock(ijk, splitFact)
m = min(ijk);
n = max(ijk) - m + 1;
p = compressPartition(partitionCartGrid(n, splitFact));

[I, J, K] = ndgrid(1 : n(1), 1 : n(2), 1 : n(3));
p         = p(ismember([I(:), J(:), K(:)],                     ...
                       ijk - m(ones([size(ijk,1), 1]), :) + 1, ...
                       'rows'));
assert(numel(p) == size(ijk,1))

%--------------------------------------------------------------------------

function ijk = getCellIJK(G)
if isfield(G.cells, 'ijkMap'),
   ijk = G.cells.ijkMap;
else
   [ijk{1:3}] = ind2sub(reshape(G.cartDims, [1, 3]), ...
                        (1 : G.cells.num) .');
   ijk = [ijk{:}];
end

%--------------------------------------------------------------------------

function s = id(s)
s = ['splitBlocks:', s];
