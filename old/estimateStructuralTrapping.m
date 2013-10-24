function volumes = estimateStructuralTrapping(G, z_spill_loc)
    vals = unique(z_spill_loc);
    vals = vals(vals>0);
    N    = numel(vals);
    volumes = zeros(N,1);
    
    for i = 1:N;
        ind = z_spill_loc == vals(i);
        volumes(i) = sum((vals(i) - G.cells.z(ind)).*G.cells.volumes(ind));
    end
end
