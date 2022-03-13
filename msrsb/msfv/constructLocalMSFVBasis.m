function basis = constructLocalMSFVBasis(CG, DG, A)
% Compute basis functions from explicit dual grid.

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
    assert(isfield(DG, 'explicit'), 'Call makeExplicitDual first!');
    
    G = CG.parent;
    nc = CG.cells.num;
    is3d = G.griddim == 3 && G.cartDims(3) > 1;
    
    n_d = numel(DG.explicit);
    
    % Preallocate the final storage 
    basis.localSol = struct('p_local', cell(nc, 1), ...
                            'cells',   cell(nc, 1));
    for i = 1:n_d
        fprintf('Solving basis for block %d of %d\n', i, n_d);
        e = DG.explicit{i};
        % Set up local and global maps for the coarse block
        inner = logicalMap(G, e.inner);
        nodes = logicalMap(G, e.nodes);
        edges = logicalMap(G, e.edges);
        faces = logicalMap(G, e.faces);
        all = inner | nodes | edges | faces;
        
        isAll = find(all);
        isInner = ismember(isAll, find(inner));
        isNodes = ismember(isAll, find(nodes));
        isEdges = ismember(isAll, find(edges));
        isFaces = ismember(isAll, find(faces));
        
        nn = sum(isNodes);
        na = numel(isAll);
        
        A_local = localSys(A, all);
        
        % Solve the local problems using successive boundary conditions
        p_n = eye(nn, nn);
        if is3d
            % 3D problems have edge -> face 
            p_e = solveLocal(A_local, isEdges, isFaces | isInner, isNodes, p_n);
            p_f = solveLocal(A_local, isFaces, isInner | isNodes, isEdges, p_e);
        else
            % 2D problems have only faces (to keep in line with mrsts
            % notion of what an interface is)
            p_f = solveLocal(A_local, isFaces, isInner, isNodes, p_n);
        end
        
        p_i = solveLocal(A_local, isInner, isEdges | isNodes, isFaces, p_f);

        % Store the pressure locally
        p = zeros(na, nn);
        p(isNodes, :) = p_n;
        p(isInner, :) = p_i;
        if is3d
            p(isEdges, :) = p_e;
        end
        p(isFaces, :) = p_f;
        p(isnan(p)) = 1;
        
        % Distribute the pressure solutions to the correct primal coarse block
        nodes = sort(e.nodes);
        for j = 1:numel(nodes)
            globalCoarse = CG.partition(nodes(j));
            basis.localSol(globalCoarse).p_local{end+1} = p(:, j);
            basis.localSol(globalCoarse).cells{end+1} = isAll;
        end
    end
    
    [i_ind, j_ind, values] = deal([]);
    for i = 1:CG.cells.num
        % Actually construct a interpolator matrix
        tmp = zeros(G.cells.num, 1);
        bloc = basis.localSol(i);
        for j = 1:numel(bloc.cells)
            % overwrite without question...
            tmp(bloc.cells{j}) = bloc.p_local{j};
        end
        subs = find(tmp ~= 0);
        
        i_ind = [i_ind; subs]; %#ok
        j_ind = [j_ind; repmat(i, numel(subs), 1)]; %#ok
        values = [values; tmp(subs)]; %#ok
    end
    basis.B = sparse(i_ind, j_ind, values, CG.parent.cells.num, CG.cells.num);
end

function p_local = solveLocal(A_local, keep, ignore, eliminate, bc)
    A = localSys(A_local, keep, ignore);
    rhs = -A_local(keep, eliminate)*bc;

    p_local = A\rhs;
end

function m = logicalMap(G, subset)
    m = false(G.cells.num, 1);
    m(subset) = true; 
end

function A_loc = localSys(A, subset, exclude)
    assert(islogical(subset));
    if nargin == 2
        exclude = ~subset;
    end
    
    if ~any(subset); A_loc = []; return; end
    
    A_loc = A(subset, subset) + diag(sum(A(subset, exclude), 2));
end
