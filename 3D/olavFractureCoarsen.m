function CG = olavFractureCoarsen(G, p_m, dof_frac, grow_dist)
    % G grid with a single fracture
    % p_m is the matrix coarse partition
    % dof_frac how many degrees of freedom in that fracture
    % grow dist is the matrix distance used to get interaction region near
    % fracture cells
    Gmc = G.Matrix.cells.num;
    
    p = zeros(G.cells.num, 1);
    p(1:Gmc) = p_m;
    pmax = max(p);
    
    % Do matrix first
    CG_m = generateCoarseGrid(G.Matrix, p_m);
    CG_m = coarsenGeometry(CG_m);
    
    CG_m = storeInteractionRegion(CG_m);
    
    % hard-coded to one fracture...
    A = getConnectivityMatrix(getNeighbourship(G.FracGrid.Frac1));
    p_f = callMetisMatrix(A, dof_frac) + pmax;
    
    p(Gmc+1 : end) = p_f;
    % Keep interaction from matrix coarse grid for those blocks, avoid
    % fractures
    interaction = cell(max(p), 1);
    for i = 1:pmax
        interaction{i} = CG_m.cells.interaction{i};
    end
    A_full = getConnectivityMatrix(getNeighbourship(G, 'Topological'));

    % Extremely hack-ish way of getting topological distance...
    for ix = 1:dof_frac
        dist = double(p == pmax + ix);
        for j = 1:grow_dist
            dist = A_full*dist;
        end
        interaction{pmax + ix} = find(dist > 0);
    end
    
    CG = generateCoarseGrid(G, p);
    CG = coarsenGeometry(CG);
    CG.cells.centers = mapCenters(CG);
    
    CG.cells.interaction = interaction;
end