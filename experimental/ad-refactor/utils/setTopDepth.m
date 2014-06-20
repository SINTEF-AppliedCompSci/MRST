function G = setTopDepth(G, depth)
% ----------------------------------------------------------------------------
    mindepth = min(G.nodes.coords(:,3));    
    reps = numel(G.nodes.coords(:,3)) ./ numel(depth);
    G.nodes.coords(:,3) = G.nodes.coords(:,3) + ...
                                  repmat((depth - mindepth), reps, 1);
end
