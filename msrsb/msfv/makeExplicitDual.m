function DG = makeExplicitDual(CG, DG)
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
    mrstModule add matlab_bgl
    require matlab_bgl
    G  = CG.parent;
    is3d = G.griddim == 3 && G.cartDims(3) > 1;
    
    N = getNeighbourship(G);
    A = getConnectivityMatrix(N);
    
    % First work out the cycles in the inner system, corresponding to
    % distinct coarse blocks
    inner = false(G.cells.num, 1);
    inner(DG.ii) = true;
    A_ii = A(inner, inner);
    c_i = components(A_ii);
    
    nc_dual = max(c_i);
    
    coarseBlocks = cell(nc_dual, 1);
    for i = 1:nc_dual
        fprintf('Processing block %d of %d...\n', i, nc_dual);
        % Then we can work out which primal coarse blocks are included in this
        % region
        local = reshape(DG.ii(c_i == i), [], 1);
        localCoarse = CG.partition(local);
        

        % Grab the correct "nodal" points
        nodes = DG.nn(unique(localCoarse));
        
        faces = neighborsByFaces(G, local, DG.ee);
        if is3d
            edges = neighborsByFaces(G, faces, DG.lineedge);
        else
            edges = [];
        end
        
        coarseBlocks{i} = struct('nodes', nodes, ...
                                'edges', edges, ...
                                'faces',  faces, ...
                                'inner', local);
    end
    DG.explicit = coarseBlocks;
end

function n = neighborsByFaces(G, subset, otherset)
    % Find neighbors of subset inside otherset
        bf = boundaryFaces(G, subset);
        n = G.faces.neighbors(bf, :);
        % Strip away boundary neighbors and place in a unique, sorted array
        n = unique(n(:));
        n = n(n~=0);
        n = n(ismember(n, otherset));
end
