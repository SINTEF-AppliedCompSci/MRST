function [D, bcvals] = setupMpsaNodeFaceBc(bc, G, tbls)
    % The structure bc gives conditions on the nodeface displacement

    assert(isfield(bc, 'bcnodefacetbl'), ['this function is meant to set ' ...
                        'boundary conditions at the nodeface'])
    
    nodefacetbl    = tbls.nodefacetbl;
    nodefacecoltbl = tbls.nodefacecoltbl;
    coltbl         = tbls.coltbl;
    facetbl        = tbls.facetbl;
    
    bcnodefacetbl = bc.bcnodefacetbl;
    bcnodefacetbl = bcnodefacetbl.addLocInd('bcinds');
    bcnodefacecoltbl = crossIndexArray(bcnodefacetbl, coltbl, {}, 'optpureproduct', true);
    
    linform = bc.linform;
    linform = reshape(linform', [], 1);
    % linform belongs to bcnodefacecoltbl    

    % We compute the (pseudo) area of the nodeface
    map = TensorMap();
    map.fromTbl = nodefacetbl;
    map.toTbl = facetbl;
    map.mergefds = {'faces'};
    map = map.setup();
    
    nnodeperface = map.eval(ones(nodefacetbl.num, 1));
    faces = facetbl.get('faces');
    nfareas = 1./nnodeperface.*(G.faces.areas(faces));
    
    % We weight the linear form with the area
    prod = TensorProd();
    prod.tbl1 = facetbl;
    prod.tbl2 = bcnodefacecoltbl;
    prod.tbl3 = bcnodefacecoltbl;
    prod.mergefds = {'faces'};
    prod = prod.setup();
    
    linform = prod.eval(nfareas, linform);
    
    prod = TensorProd();
    prod.tbl1 = facetbl;
    prod.tbl2 = bcnodefacetbl;
    prod.tbl3 = bcnodefacetbl;
    prod.mergefds = {'faces'};
    prod = prod.setup();
    
    bcvals = bc.linformvals;
    bcvals = prod.eval(nfareas, bcvals);
    
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


