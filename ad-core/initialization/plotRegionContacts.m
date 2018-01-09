function plotRegionContacts(G, region)
    coord = G.nodes.coords;
    xmin = min(coord(:, 1));
    xmax = max(coord(:, 1));
    ymin = min(coord(:, 2));
    ymax = max(coord(:, 2));
    X = [xmax, xmin, xmin, xmax];
    Y = [ymax, ymax, ymin, ymin];
    nc = numel(region.contacts);
    colors = jet(nc);


    for i = 1:nc
        c = region.contacts(i);
        Z = [c, c, c, c];
        patch(X, ...
              Y, ...
              Z, colors(i, :));
    end
end