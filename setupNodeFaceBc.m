function [D, bcvals] = setupNodeFaceBc(bc, G, tbls)
% the structur bc gives conditions on the nodeface displacement

    assert(isfield(bc, 'bcnodefacetbl'), ['this function is meant to set ' ...
                        'boundary conditions at the nodeface'])
    
    nodefacetbl    = tbls.nodefacetbl;
    nodefacecoltbl = tbls.nodefacecoltbl;
    coltbl         = tbls.coltbl;
    
    bcnodefacetbl = bc.bcnodefacetbl;
    bcnodefacetbl = bcnodefacetbl.addLocInd('bcinds');
    bcnodefacecoltbl = crossTable(bcnodefacetbl, coltbl, {}, 'optpureproduct', true);
    
    linform = bc.linform;
    linform = reshape(linform', [], 1);
    % linform belongs to bcnodefacecoltbl    
    
    bcvals = bc.linformvals;
    
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


