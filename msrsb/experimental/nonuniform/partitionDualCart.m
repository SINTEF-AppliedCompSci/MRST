function CG = partitionDualCart(CG)
%Undocumented Utility Function

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    if G.griddim == 3 && G.cartDims(3) > 1
        dimension = 3;
    else
        dimension = 2;
    end
    
    ijk = gridLogicalIndices(G);
    
    counts = zeros(G.cells.num, 1);
    centers = CG.cells.centers;
    for i = 1:dimension
        li = ijk{i};
        
        tmp = false(max(li), 1);
        tmp(li(centers)) = true;
        
        counts = counts + tmp(li);
    end
    
    if dimension == 2
        DG.nn = centers;
        DG.lineedge = [];
        DG.ee = find(counts == 1);
        DG.ii = find(counts == 0);
    else 
        DG.nn = centers;
        DG.lineedge = find(counts == 2);
        DG.ee = find(counts == 1);
        DG.ii = find(counts == 0);
    end
    CG.dual = makeExplicitDual(CG, DG);
end
