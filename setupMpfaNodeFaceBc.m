function [D, bcvals] = setupMpfaNodeFaceBc(bc, tbls)

    bcnodefacetbl = bc.bcnodefacetbl;
    bcvals        = bc.bcvals;
    nodefacetbl   = tbls.nodefacetbl;
    
    map = TensorMap();
    map.fromTbl = bcnodefacetbl;
    map.toTbl = nodefacetbl;
    map.mergefds = {'nodes', 'faces'};
    map = map.setup();
    
    D_T = SparseTensor();
    D_T = D_T.setFromTensorMap(map);
    D = D_T.getMatrix();
    
    
end


