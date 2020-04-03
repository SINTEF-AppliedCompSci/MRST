function op = setupCell2NodeDispOperator(G, tbls)
    
    nodetbl = tbls.nodetbl;
    celltbl = tbls.celltbl;
    cellcoltbl = tbls.cellcoltbl;
    nodecoltbl = tbls.nodecoltbl;
    cellnodetbl = tbls.cellnodetbl;
    cellnodecoltbl = tbls.cellnodecoltbl;
    
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = nodetbl;
    map.mergefds = {'nodes'};
    map = map.setup();
    
    coef = map.eval(ones(cellnodetbl.num, 1));
    coef = 1./coef;
    
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
    
    op = SparseTensor();
    op = op.setFromTensorProd(coef, prod);
    
    op = op.getMatrix();
    
end
