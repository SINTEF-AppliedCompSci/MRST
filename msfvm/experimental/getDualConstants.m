function DG = getDualConstants(CG, DG)
%Undocumented Utility Function

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

    G = CG.parent;
    is3d = CG.griddim == 3 && G.cartDims(3) > 1;
    [~, DG] = createPermutationMatrix(sparse(G.cells.num, G.cells.num), DG, CG, 'Speedup', is3d);

    P = double(DG.P);
    X = restrictOperator(CG, DG, DG.N);
    DG.X = X*P;
end

function X = restrictOperator(CG, DG, Nf)
    %Returns a matrix representing the restriction operator
    %Should be nxf big where n is coarse nodes and f total fine nodes
    %Permute the partition to the new index space
    permuted_partition = DG.P*CG.partition;
    xind = zeros(1,Nf);
    yind = zeros(1,Nf);
    pos = 1;
    for coarse = 1:CG.cells.num
       %find indices of corresponding coarse nodes
       indices = find(permuted_partition == coarse);
       %insert the correct cells into the arrays of indices
       M = length(indices);
       yind(pos:(pos + M-1)) = indices;
       xind(pos:(pos + M-1)) = coarse;
       pos = pos + M;
    end
    X = sparse(xind, yind, 1, Nf, Nf) > 0;
    return
end
