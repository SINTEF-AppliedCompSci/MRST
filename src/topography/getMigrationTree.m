function subtraps = getMigrationTree(G, A, trap, depth)
% Recursively traverse and find the full migration tree
    subtraps = find(A(trap, :));
    tmp = [];
    for i = 1:numel(subtraps);
        trp = getMigrationTree(G, A, subtraps(i), depth + 1);
        tmp = [tmp trp]; %#ok
    end
    subtraps = [reshape(subtraps, 1, []) reshape(tmp, 1, [])];
end
