function cont = extractGridBoundaryNodes(Gt)
%  NB: It is assumed that there are no internal holes, i.e. that there is one
%  unique boundary.

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

    % Algorithm depends on each edge having exactly two nodes
    assert(unique(diff(Gt.faces.nodePos)) == 2);
    
    % Identify end nodes of all edges
    edgenodes = reshape(Gt.faces.nodes, 2, [])';
        
    % Determine indices of boundary edges
    bnd_edges = find(prod(Gt.faces.neighbors,2)==0);
    
    % Identify endnodes of boundary edges, and sort for each edge
    bnd_nodes = sort(edgenodes(bnd_edges,:), 2, 'ascend');%#ok
    
    max_node_ix = max(bnd_nodes(:));
    
    % Establishing upper triangular part of adjacency matrix
    A = sparse(bnd_nodes(:,1), bnd_nodes(:,2), 1, max_node_ix, max_node_ix);
    
    % Assembling result (@ brute force loop - perhaps a more elegant solution exists?)

    cont = max_node_ix; % We choose this as the first node of our loop
    next = 0;
    while next ~= cont(1)
        next = find(A(:,cont(end)), 1); 
        if isempty(next)
            next = find(A(cont(end), :), 1);
            A(cont(end), next) = 0; %#ok
        else
            A(next, cont(end)) = 0; %#ok
        end
        cont = [cont; next]; %#ok
    end
end
