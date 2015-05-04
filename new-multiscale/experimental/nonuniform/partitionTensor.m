function p = partitionTensor(G, di, dj, dk)
    if G.griddim == 3
        delta = {di, dj, dk};
        assert(nargin == 4);
    else
        delta = {di, dj};
        assert(nargin == 3);
    end
    ijk = gridLogicalIndices(G);
    ijk = [ ijk{:} ];
    M = max(ijk) - min(ijk) + 1;
    
    assert(all(cellfun(@sum, delta) == M));
    
    
    
    p = ones(G.cells.num, 1);
    for i = G.griddim:-1:1
        d = delta{i};
        
        p_loc = rldecode((1:numel(d))', d);
        
        p = p + p_loc(ijk(:, i) - min(ijk(:, i)) + 1)*max(p);
        p = processPartition(G, p);
        p = compressPartition(p);
    end
end
