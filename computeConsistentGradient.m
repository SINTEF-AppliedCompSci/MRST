function g = computeConsistentGradient(G, eta, tbls, mappings, varargin)
    
    opt = struct('bcetazero', false);
    opt = merge_options(opt, varargin{:});
    bcetazero = opt.bcetazero;
    
    cellnodefacecents = computeNodeFaceCentroids(G, tbls, eta, 'bcetazero', opt.bcetazero);

    coltbl             = tbls.coltbl;
    cellnodetbl        = tbls.cellnodetbl;
    cellnodecoltbl     = tbls.cellnodecoltbl;
    cellnodefacetbl    = tbls.cellnodefacetbl;
    cellnodefacecoltbl = tbls.cellnodefacecoltbl;
    
    cellnode_from_cellnodeface = mappings.cellnode_from_cellnodeface;
    
    d_num    = coltbl.num;
    cn_num   = cellnodetbl.num;
    cnc_num  = cellnodecoltbl.num; 
    cnfc_num = cellnodefacecoltbl.num;
    cnf_num  = cellnodefacetbl.num;
    
    [c, i] = ind2sub([d_num, cnf_num], (1 : cnfc_num)');
    ind1 = i;
    ind2 = sub2ind([d_num, cn_num], c, cellnode_from_cellnodeface(i));

    assert(cnc_num == cnf_num, ['This implementation of mpsaw cannot handle ' ...
                        'this grid']);

    A = sparse(ind1, ind2, cellnodefacecents, cnc_num, cnc_num);

    opt.invertBlocks = 'mex';
    bi = blockInverter(opt);

    sz = repmat(coltbl.num, cellnodetbl.num, 1);
    invA = bi(A, sz);

    ind = sub2ind([cnf_num, cnc_num], ind2, ind1);
    
    g = invA(ind);
    
end

