function CG = basisToCoarseGrid(G, I)
% Create a coarse grid from a set of basis functions so that each coarse
% block corresponds to the places where that basis is the largest value.
    for i = 1:size(I, 2)
        [~, biggest] = max(I(:, i));
        I(biggest, i) = inf;
    end
    [~, p] = max(I, [], 2);
    
    CG = generateCoarseGrid(G, p);
    CG = coarsenGeometry(CG);
end