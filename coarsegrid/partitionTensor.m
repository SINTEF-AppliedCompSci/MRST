function p = partitionTensor(G, di, dj, dk)
%Partition Logically Cartesian Grid Into Tensor Product Blocks
%
% SYNOPSIS:
%   % In two space dimensions
%   p = partitionTensor(G, di, dj)
%
%   % In three space dimensions
%   p = partitionTensor(G, di, dj, dk)
%
% PARAMETERS:
%   G          - MRST grid as described in grid_structure.  Assumed to
%                carry an underlying logically Cartesian structure (i.e.,
%                each cell can be identified by its (I,J) pair or (I,J,K)
%                triplet in two and three space dimensions, respectively.
%
%   di,dj,dk   - Tensor product block sizes.  Vectors of positive integers.
%                Specifically, the block at position (x,y,z) consists of
%
%                      di(x)-by-dj(y)-by-dk(z)
%
%                cells from the underlying fine-scale grid (i.e., G)--at
%                least if all cells are active.  These inputs play the same
%                roles, and are subject to the same restrictions like
%
%                      ALL([SUM(di), SUM(dj), SUM(dk)] == G.cartDims)
%
%                as the 'dimDist' inputs of function MAT2CELL.
%
% RETURNS:
%   p - Partition vector.  Contains no empty or multiply connected blocks.
%
% SEE ALSO:
%   `partitionUI`, `partitionLayers`, `mat2cell`.

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

    if G.griddim == 3
        assert(nargin == 4, ...
              ['Function must be called with three block-size ', ...
               'inputs in 3D']);
        delta = {di, dj, dk};
    else
        assert(nargin == 3, ...
              ['Function must be called with two block-size ', ...
               'inputs in 2D']);
        delta = {di, dj};
    end
    ijk = gridLogicalIndices(G);
    ijk = [ ijk{:} ];
    M = max(ijk) - min(ijk) + 1;

    assert(all(cellfun(@sum, delta) == M));



    p = ones(G.cells.num, 1);
    for i = G.griddim:-1:1
        d = delta{i};

        p_loc = rldecode((1:numel(d))', d);

        p = p + p_loc(ijk(:, i) - min(ijk(:, i)) + 1)*max(p);
        p = processPartition(G, p);
        p = compressPartition(p);
    end
end
