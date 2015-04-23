function n = getCellNeighbors(G, c)
    n = G.faces.neighbors(any(G.faces.neighbors == c, 2), :);
    n = n(:);
    n = n(n~=0);
    n = n(n~=c);
    n = unique(n);
end