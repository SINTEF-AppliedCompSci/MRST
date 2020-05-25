function faces = addSealingFaces(G, varargin)
    opt = struct('i_range', [-inf, inf], ...
                 'j_range', [-inf, inf], ...
                 'k_range', [-inf, inf], ...
                 'x_range', [-inf, inf], ...
                 'y_range', [-inf, inf], ...
                 'z_range', [-inf, inf], ...
                 'fraction', 1 ...
                 );

    opt = merge_options(opt, varargin{:});
    [ii, jj, kk] = gridLogicalIndices(G);
    neighbors = G.faces.neighbors;
    act = all(neighbors > 0, 2);
    N = neighbors(act, :);
    all_faces = find(act);
    
    
    x = G.faces.centroids(act, 1);
    y = G.faces.centroids(act, 2);
    z = G.faces.centroids(act, 3);
    
    log_mask = all(ii(N) >= opt.i_range(1), 2) & all(ii(N) <= opt.i_range(2), 2) & ...
               all(jj(N) >= opt.j_range(1), 2) & all(jj(N) <= opt.j_range(2), 2) & ...
               all(kk(N) >= opt.k_range(1), 2) & all(kk(N) <= opt.k_range(2), 2);
    
    coord_mask =  all(x >= opt.x_range(1), 2) & all(x <= opt.x_range(2), 2) & ...
                  all(y >= opt.y_range(1), 2) & all(y <= opt.y_range(2), 2) & ...
                  all(z >= opt.z_range(1), 2) & all(z <= opt.z_range(2), 2);
    mask = log_mask & coord_mask;
    mask = mask & kk(N(:, 1)) ~= kk(N(:, 2));
    
    faces = all_faces(mask);
    if opt.fraction < 1
        tmp = rand(size(faces));
        faces = faces(tmp < opt.fraction);
    end
end

