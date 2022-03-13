function M = subregionConnectivityMatrix(Gt, nodes, use_diags)
%
% Generate the connectivity matrix between a subset of the nodes of a 2D grid.
%
% SYNOPSIS:
%   function M = subregionConnectivityMatrix(Gt, nodes)
%
% DESCRIPTION:
%
% PARAMETERS:
%   Gt        - top surface grid
%   nodes     - vector with node indices, representing a subset of the nodes in
%               Gt.
%   use_diags - if 'true', cell diagonals will also count as connections
% 
% RETURNS:
%   M - Connectivity matrix between the nodes in 'nodes', as defined by the
%       edges of the top surface grid Gt, or if diagonals are used, all nodes
%       sharing a cell are considered connected.
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
  
   if use_diags
      neigh_rels = cell_node_neighs(Gt, nodes);
   else
      neigh_rels = edge_node_neighs(Gt, nodes);
   end
      
    dists = compute_dists(Gt, neigh_rels);
    
    % establish connectivity matrix
    n = numel(nodes);
    reindex = zeros(max(nodes),1);
    reindex(nodes) = 1:numel(nodes);
    M = sparse(reindex(neigh_rels(:,1)), ...
               reindex(neigh_rels(:,2)), ...
               dists, n, n);
    M = M+M';
end

function nrels = edge_node_neighs(Gt, nodes)
    % list of all node connectivities
    neigh_rels = reshape(Gt.faces.nodes, 2, [])';
    
    % elminiate connectivities not related to the specified nodes
    ixs = ismember(neigh_rels(:,1), nodes);
    neigh_rels = neigh_rels(ixs,:);
    ixs = ismember(neigh_rels(:,2), nodes);
    nrels = neigh_rels(ixs,:);
end


function d = compute_dists(Gt, neigh_rels)

    c1 = Gt.nodes.coords(neigh_rels(:,1),:);
    c2 = Gt.nodes.coords(neigh_rels(:,2),:);
    
    d = sqrt(sum((c1 - c2).^2, 2));
end

function m = cell_nodes_map(Gt, nodes)

   cn = getSortedCellNodes(Gt);
   icells = rldecode((1:Gt.cells.num)', diff(Gt.cells.facePos));
   ix = ismember(cn, nodes);
   
   m = sparse(icells(ix), cn(ix), 1, Gt.cells.num, Gt.nodes.num);
   
end

function nrels = cell_node_neighs(Gt, nodes)
   
   m = cell_nodes_map(Gt, nodes);
   ncells = find(sum(m,2) > 1); % cells with more than 1 node from the set
   nrels = arrayfun(@(c) construct_clique(find(m(c, :))), ncells, 'uniformoutput', false);
   nrels = vertcat(nrels{:});
end

function res = construct_clique(nodes)

   res = nodes(nchoosek(1:numel(nodes), 2));
end
