function D = setupBC(bcstruct, G, tbls)

    
    n = numel(bcstruct);
    D = cell(n, 1);
    
    for i = 1 : n
        extfaces = bcstruct{i}.extfaces;
        linform  = bcstruct{i}.linform;
        D{i} = assignLinearForm(extfaces, linform, tbls);
    end
    
    D = horzcat(D{:});
    
end

function D = assignLinearForm(extfaces, linform, tbls)
    
    nodefacetbl    = tbls.nodefacetbl;
    nodefacecoltbl = tbls.nodefacecoltbl;
    coltbl         = tbls.coltbl;
  
    extfacetbl.faces = extfaces;
    extfacetbl = IndexTable(extfacetbl);

    extnodefacetbl = crossTable(nodefacetbl, extfacetbl, {'faces'});

    prod = TensorProd();
    prod.tbl1 = coltbl;
    prod.tbl2 = extnodefacetbl;
    prod.tbl3 = nodefacecoltbl;
    prod = prod.setup();

    D_T = SparseTensor('matlabsparse', true);
    D_T = D_T.setFromTensorProd(linform, prod);
    D = D_T.getMatrix();
    
end
