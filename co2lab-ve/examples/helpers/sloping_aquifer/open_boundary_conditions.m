function bc = open_boundary_conditions(G, pfun, cross_sectional_problem)

    bc = []; % Start with an empty set of boundary faces
    
    if G.griddim == 3
        vface_ind = (G.faces.normals(:,3) == 0); % identify all vertical faces
        depth = G.faces.centroids(:, 3);
    else
        assert(G.griddim == 2); % presumably a VE grid
        vface_ind = true(size(G.faces.normals, 1), 1); % all faces are vertical
        depth = G.faces.z;
    end
    
    bface_ind = (prod(G.faces.neighbors, 2) == 0); % identify all boundary faces 
    bc_face_ix = find(vface_ind & bface_ind); % identify all lateral boundary faces

    if cross_sectional_problem
        remove = G.faces.normals(bc_face_ix, 2) ~= 0 | ...
                 G.faces.centroids(bc_face_ix, 1) == 0;
        bc_face_ix(remove) = [];
    end
    
    
    
    p_face_pressure = pfun(depth(bc_face_ix));
    
    % p_face_pressure = cell_pressures(bc_cell_ix);  % set lateral boundary face
    %                                                % pressure equal to pressure
    %                                                % of corresponding cell

    bc = addBC(bc, ...
               bc_face_ix, 'pressure', ...      % Add hydrostatic pressure         
               p_face_pressure, 'sat', [1, 0]); % conditions to open boundary faces
end
