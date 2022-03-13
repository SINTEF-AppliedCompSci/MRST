function blockIx = partitionLayers(G, coarseDim, L)
%Partition grid uniformly in logical (I,J) direction, non-uniformly in K.
%
% SYNOPSIS:
%   blockIx = partitionLayers(G, coarseDim, L)
%
% DESCRIPTION:
%   Partitions a corner point grid in grid_structure format (relatively)
%   non-uniformly in k-direction and uniformly in ij-directions in index
%   (logical) space.  In the ij-directions, each coarse block is
%   comprised of roughly the same number of cells from the fine grid.
%
% PARAMETERS:
%   G         - grid_structure structure having either a valid field
%               'G.cells.ijkMap' or both valid fields 'G.cartDims' and
%               'G.cells.indexMap'.
%
%   coarseDim - Number of coarse blocks in each ij direction.
%               Assumed to be a LENGTH 2 vector of integers.
%
%   L         - A run-length-encoded vector of the desired layers.
%
% RETURNS:
%   blockIx   - A G.cells.num-by-1 vector mapping cells to coarse blocks.
%
% EXAMPLE:
%   See simpleCornerPointExampleMS
%
% SEE ALSO:
%   `processPartition`, `partitionUI`.

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


if ~grid_ok(G)
   error('partitionLayers:grid_structure:NOIJK', ...
         'Grid is not a valid grid_structure structure');
end

if numel(coarseDim) ~= 2
   error('partitionLayers:coarseDim:Invalid', ...
         '`coarseDim'' is not valid coarse ij dimension spec');
end

if isfield(G.cells, 'ijkMap')
   ijk = double(G.cells.ijkMap);
   M   = double(max(ijk));
else
   [ijk{1:3}] = ind2sub(G.cartDims, G.cells.indexMap);
   ijk = double([ijk{:}]);
   M   = double(G.cartDims);
end

%blockIx = zeros([G.cells.num, 1]);
L    = reshape(L, 1, []);
lNum = rldecode(1 : numel(L) - 1, diff(L), 2) - 1;
blockIx = lNum(ijk(:,3))';
for d = 2 : -1 : 1
   B = coarseDim(d);
   blockIx = lbLinDist(ijk(:,d) - 1, M(d), B) + B*blockIx;
end
blockIx = blockIx + 1;  % Map block numbers to (1..PROD(coarseDim))


function f = lbLinDist(f, M, B)
% lbLinDist -- Load-balanced linear distribution
% See Eric F. Van de Velde, Concurrent Scientific Computing,
% 1994, Springer Verlag, p. 54 (Sect. 2.3) for details.
%
% Maps index set (0..M-1) to blocks (0..B-1).

L = floor(M ./ B);  % Tentative number of cells per coarse block.
R = mod(M, B);      % Additional cells not previously accounted for.
f = max(floor(f ./ (L + 1)), floor((f - R) ./ L));


function bool = grid_ok(G)
bool = ~isempty(G)          && isstruct(G)       && ...
        isfield(G, 'cells') &&                      ...
       ~isempty(G.cells)    && isstruct(G.cells);
bool = bool && ((isfield(G.cells, 'ijkMap') && ...
                 size(G.cells.ijkMap,1) == G.cells.num) || ...
                (isfield(G, 'cartDims') && isfield(G.cells, 'indexMap')));
