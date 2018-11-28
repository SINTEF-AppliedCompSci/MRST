function GH = refineGrid(GC, GF, cells)

    GH = GC;
    ix = any(GH.partition == cells',2);
    partition     = GC.partition;
    partition     = partition + nnz(ix) - 1;
    partition(ix) = 1:nnz(ix);
    
    GH = generateCoarseGrid(GF, partition);
    GH = coarsenGeometry(GH);
    
end