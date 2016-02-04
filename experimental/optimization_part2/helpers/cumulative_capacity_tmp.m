function cumcap = cumulative_capacity_tmp(trapcap, adj)
% taken from optimizeFormation.m

    % Compute cumulative capacity of each trap and those it is connected to upstream
    cumcap = trapcap;
    adj_orig = adj;
    while nnz(adj)
        cumcap = cumcap + adj * trapcap;
        adj = spones(adj * adj_orig);
    end
end


