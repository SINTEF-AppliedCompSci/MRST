function blockIx = partitionUI(G, coarseDim)
%Partition grid uniformly in logical space.
%
% SYNOPSIS:
%   blockIx = partitionUI(G, coarseDim)
%
% DESCRIPTION:
%   Partitions a corner point grid in grid_structure format (relatively)
%   uniformly in index (logical) space.  Each coarse block is
%   comprised of roughly the same number of cells from the fine grid.
%
% PARAMETERS:
%   G         - grid_structure structure having valid 'G.cartDims' and
%               'G.cells.indexMap' fields.
%
%   coarseDim - Number of coarse blocks in each physical direction.
%               Assumed to be a LENGTH 2 or 3 vector of integers.
%
% RETURNS:
%   blockIx   - A G.cells.num-by-1 vector mapping cells to coarse blocks.
%
% SEE ALSO:
%   `processPartition`, `partitionLayers`.

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

if ~isfield(G, 'cartDims')
   error(msgid('grid_structure:Invalid'), ...
         ['Grid must be be logically Cartiesian to use partitionUI.'...
         ' Missing field ''cartDims''.']);
end
if ~grid_ok(G)
   error(msgid('grid_structure:Invalid'), ...
         'Grid is not a valid grid_structure structure');
end

if numel(coarseDim) ~= G.griddim || any(coarseDim < 1)
   error(msgid('coarseDim:Invalid'), ...
        ['Parameter ''coarseDim'' must be a %d component ', ...
         'vector\ncontaining strictly positive integers.'], G.griddim);
end

assert(all(coarseDim > 0), ...
   'Negative element(s) in parameter ''coarseDim''.');
assert(all(coarseDim <= G.cartDims), ...
   'Parameter ''coarseDim'' exceed G.cartDims.');

[i{1 : numel(G.cartDims)}] = ...
   ind2sub(double(G.cartDims), double(G.cells.indexMap));
i = [ i{:} ];
M = max(i, [], 1) - min(i, [], 1) + 1;

blockIx = zeros([G.cells.num, 1]);
for d = numel(coarseDim) : -1 : 1
   B = double(coarseDim(d));
   blockIx = lbLinDist(i(:,d) - min(i(:,d)), M(d), B) + B*blockIx;
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
bool = bool && isfield(G, 'cartDims') && isfield(G.cells, 'indexMap');
