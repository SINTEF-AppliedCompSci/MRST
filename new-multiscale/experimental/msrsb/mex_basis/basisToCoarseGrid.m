function CG = basisToCoarseGrid(G, I)
    for i = 1:size(I, 2)
        [~, biggest] = max(I(:, i));
        I(biggest, i) = inf;
    end
    [~, p] = max(I, [], 2);
    
    CG = generateCoarseGrid(G, p);
    CG = coarsenGeometry(CG);
end