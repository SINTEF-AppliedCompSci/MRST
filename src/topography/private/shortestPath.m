function path = shortestPath(M, start, target)

    num_nodes = size(M,1);
    assert(size(M,2) == num_nodes);
    path = target;
    
    % Initialize distance vector and neighbor vector
    D = ones(num_nodes, 1) * Inf;
    n = zeros(num_nodes,1);
    D(start) = 0;
    
    % @@ Currently only intended to work on relatively small datasets, so we
    % allow ourselves to work with the full representation of M
    M = full(M);
    M(M==0) = inf;
    
    while isinf(D(target))

        finished_ixs   = find(~isinf(D));
        unfinished_ixs = find(isinf(D));

        tmp = M(finished_ixs,:);
        tmp = tmp(:, unfinished_ixs);
        
        if all(isinf(tmp(:)))
            % no pathway exits.  the nodes 'start' and 'target' are not
            % connected
            return
        end
        
        tmp = bsxfun(@plus, D(finished_ixs), tmp);
        
        [val, ix] = min(tmp,[], 1);
        
        found = isfinite(val);
        assert(~all(~found)); % should be a consequence of the previous 'if' clause
        
        D(unfinished_ixs(found)) = val(found);
        n(unfinished_ixs(found)) = finished_ixs(ix(found));
    end
    
    while path(end) ~= start
        path = [path; n(path(end))];
    end
    
    path = flipud(path);
end
