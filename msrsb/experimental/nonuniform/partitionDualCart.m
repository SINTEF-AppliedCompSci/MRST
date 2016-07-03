function CG = partitionDualCart(CG)

    G = CG.parent;

    if G.griddim == 3 && G.cartDims(3) > 1
        dimension = 3;
    else
        dimension = 2;
    end
    
    ijk = gridLogicalIndices(G);
    
    counts = zeros(G.cells.num, 1);
    centers = CG.cells.centers;
    for i = 1:dimension
        li = ijk{i};
        
        tmp = false(max(li), 1);
        tmp(li(centers)) = true;
        
        counts = counts + tmp(li);
    end
    
    if dimension == 2
        DG.nn = centers;
        DG.lineedge = [];
        DG.ee = find(counts == 1);
        DG.ii = find(counts == 0);
    else 
        DG.nn = centers;
        DG.lineedge = find(counts == 2);
        DG.ee = find(counts == 1);
        DG.ii = find(counts == 0);
    end
    CG.dual = makeExplicitDual(CG, DG);
end