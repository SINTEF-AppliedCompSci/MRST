function D = setupBC(bc, G, tbls)
    
    nodefacetbl    = tbls.nodefacetbl;
    nodefacecoltbl = tbls.nodefacecoltbl;
    coltbl         = tbls.coltbl;

    bcfacetbl.faces = bc.extfaces;
    bcfacetbl = IndexTable(bcfacetbl);
    bcfacetbl = bcfacetbl.addLocInd('bcinds');
    bcfacecoltbl = crossTable(bcfacetbl, coltbl, {}, 'optpureproduct', true);
    
    linform = bc.linform;
    linform = reshape(linform', [], 1);
    % linform belongs to bcfacecoltbl    
    
    bcnodefacetbl = crossTable(bcfacetbl, nodefacetbl, {'faces'});
    bcnodefacecoltbl = crossTable(bcnodefacetbl, coltbl, {}, 'optpureproduct', true);
    
    map = TensorMap();
    map.fromTbl = bcfacecoltbl;
    map.toTbl = bcnodefacecoltbl;
    map.mergefds = {'faces', 'bcinds', 'coldim'};
    map = map.setup();
    
    linform = map.eval(linform);
    % linform now belongs to bcnodefacecoltbl
    
    prod = TensorProd();
    prod.tbl1 = bcnodefacecoltbl;
    prod.tbl2 = bcnodefacecoltbl;
    prod.tbl3 = bcnodefacetbl;
    prod.mergefds = {'nodes', 'faces', 'bcinds'};
    prod.reducefds = {'coldim'};
    prod = prod.setup();
     
    D_T = SparseTensor('matlabsparse', true);
    D_T = D_T.setFromTensorProd(linform, prod);

    map = TensorMap();
    map.fromTbl = nodefacecoltbl;
    map.toTbl = bcnodefacecoltbl;
    map.mergefds = {'nodes', 'faces', 'coldim'};
    map = map.setup();
    
    M_T = SparseTensor('matlabsparse', true);
    M_T = M_T.setFromTensorMap(map);
    
    D_T = D_T*M_T;
    
    D = D_T.getMatrix();
    D = D';
    
end


