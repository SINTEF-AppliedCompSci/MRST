function CG = storeInteractionRegionCart(CG)
    % For when grid was created by partitionUI
    G = CG.parent;
    
    [cc, fc] = mapCenters(CG);
    if ~isfield(CG.cells, 'centers')
        CG.cells.centers = cc;
    end
    
    if ~isfield(CG.faces, 'centers')
        CG.faces.centers = fc;
    end
    
    ijk = gridLogicalIndices(G);
    ijk = [ijk{:}];
    
    interaction = cell(CG.cells.num, 1);
    for i = 1:CG.cells.num
        n = coarseNeighbors(CG, i, false);
        ijk_near = ijk(CG.cells.centers(n), :);
        
        ijkloc = ijk(CG.partition == i, :);
        
        M_ix = max(ijk_near, [], 1);
        m_ix = min(ijk_near, [], 1);
        
        M_ix = max(M_ix, max(ijkloc + 1, [], 1));
        m_ix = min(m_ix, min(ijkloc - 1, [], 1));
        
        ok = true(G.cells.num, 1);
        for j = 1:G.griddim
            ok = ok & ijk(:, j) < M_ix(j);
            ok = ok & ijk(:, j) >  m_ix(j);
        end
        interaction{i} = find(ok);
    end
    CG.cells.interaction = interaction;
end
