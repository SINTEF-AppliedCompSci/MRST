function M = subregionConnectivityMatrix(Gt, nodes)

    % list of all node connectivities
    neigh_rels = reshape(Gt.faces.nodes, 2, [])';
    
    % elminiate connectivities not related to the specified nodes
    ixs = ismember(neigh_rels(:,1), nodes);
    neigh_rels = neigh_rels(ixs,:);
    ixs = ismember(neigh_rels(:,2), nodes);
    neigh_rels = neigh_rels(ixs,:);
    
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


function d = compute_dists(Gt, neigh_rels)

    c1 = Gt.nodes.coords(neigh_rels(:,1),:);
    c2 = Gt.nodes.coords(neigh_rels(:,2),:);
    
    d = sqrt(sum((c1 - c2).^2, 2));
    
end
