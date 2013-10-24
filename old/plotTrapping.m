function plotTrapping(G, z_spill_loc, volumes)
    oldmap = colormap();
    vals = unique(z_spill_loc);
    vals = vals(vals>0);
    N = numel(vals);
    map = colormap(hsv(N));
    
    [v ind] = sort(volumes);
    for i = 1:N
        j = ind(i);
        plotGrid(G, z_spill_loc == vals(j), 'facecolor', map(i,:), 'edgec', 'none')
    end
    colormap(oldmap);
end
