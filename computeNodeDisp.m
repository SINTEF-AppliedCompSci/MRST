function unode = computeNodeDisp(ucell, tbls)
    celltbl        = tbls.celltbl;
    cellcoltbl     = tbls.cellcoltbl;
    nodetbl        = tbls.nodetbl;
    nodecoltbl     = tbls.nodecoltbl;
    cellnodetbl    = tbls.cellnodetbl;
    cellnodecoltbl = tbls.cellnodecoltbl;
    
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = nodetbl;
    map.mergefds = {'nodes'};
    map = map.setup();
    
    ncpernode = map.eval(ones(cellnodetbl.num, 1));
    coef = 1./ncpernode;
    
    map = TensorMap();
    map.fromTbl = nodetbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'nodes'};
    map = map.setup();
    
    coef = map.eval(coef);
    
    prod = TensorProd();
    prod.tbl1 = cellnodetbl;
    prod.tbl2 = cellcoltbl;
    prod.tbl3 = nodecoltbl;
    prod.reducefds = {'cells'};
    prod = prod.setup();
    
    unode = prod.eval(coef, ucell);

end

