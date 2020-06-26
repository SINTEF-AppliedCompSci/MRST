function aver = cellAverageOperator(tbls, mappings)
    
    cellnodecolrowtbl = tbls.cellnodecolrowtbl;
    cellcolrowtbl = tbls.cellcolrowtbl;
    cellnodetbl = tbls.cellnodetbl;
    celltbl = tbls.celltbl;
    coltbl = tbls.coltbl;
    
    cell_from_cellnode = mappings.cell_from_cellnode;
    
    % shortcuts
    d_num = coltbl.num;
    c_num = celltbl.num;
    cn_num = cellnodetbl.num;
    cncr_num = cellnodecolrowtbl.num;
    
    % Compute number of nodes per cells
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = celltbl;
    map.mergefds = {'cells'};
    map = map.setup();
    
    nnodepercell = map.eval(ones(cellnodetbl.num, 1));
    
    % Compute cell average stress
    prod = TensorProd();
    prod.tbl1 = celltbl;
    prod.tbl2 = cellnodecolrowtbl;
    prod.tbl3 = cellcolrowtbl;
    prod.mergefds = {'cells'};
    
    prod.pivottbl = cellnodecolrowtbl;
    [r, c, i] = ind2sub([d_num, d_num, cn_num], (1 : cncr_num)');
    prod.dispind1 = cell_from_cellnode(i);
    prod.dispind2 = (1 : cncr_num)';
    prod.dispind3 = sub2ind([d_num, d_num, c_num], r, c, cell_from_cellnode(i));
    prod.issetup = true;

    aver_T = SparseTensor();
    aver_T = aver_T.setFromTensorProd(1./nnodepercell, prod);

    aver = aver_T.getMatrix();
    
end