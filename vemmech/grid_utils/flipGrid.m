function G_new = flipGrid(G);
    G_new = G;
    per = [3, 1, 2];
    tags = reshape(1:6, 2, []);
    new_tag = reshape(tags(:, [2, 3, 1]), [], 1);
    G_new.nodes.coords = G_new.nodes.coords(:, per);
    G_new.cells.faces(:, 2) = new_tag(G.cells.faces(:, 2));
end