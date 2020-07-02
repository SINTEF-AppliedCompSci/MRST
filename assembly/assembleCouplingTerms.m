function assembly = assembleCouplingTerms(G, eta, alpha, nnodespercell, tbls, mappings)
    
    cellnodefacecoltbl = tbls.cellnodefacecoltbl;
    nodefacecoltbl     = tbls.nodefacecoltbl;
    cellnodefacetbl    = tbls.cellnodefacetbl;
    cellnodetbl        = tbls.cellnodetbl;
    cellcoltbl         = tbls.cellcoltbl;
    celltbl            = tbls.celltbl;
    
    % We fetch the vector g, which belongs to cellnodefacecoltbl and is used to
    % construct the consistent divergence operator.
    g = computeConsistentGradient(G, eta, tbls, mappings);

    % We fetch the vector facetNormals, which belongs to cellnodefacecoltbl and is
    % used to construct the finite volume divergence operator.
    normals = computeFacetNormals(G, cellnodefacetbl);

    % Multiply with Biot's coefficient alpha
    prod = TensorProd();
    prod.tbl1 = celltbl;
    prod.tbl2 = cellnodefacecoltbl;
    prod.tbl3 = cellnodefacecoltbl;
    prod.mergefds = {'cells'};
    prod = prod.setup();
    
    normals = prod.eval(alpha, normals);
    
    %% We setup the finite volume divergence operator
    % divfv : nodefacecoltbl -> celltbl
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacecoltbl;
    prod.tbl2 = nodefacecoltbl;
    prod.tbl3 = celltbl;
    prod.reducefds = {'nodes', 'faces', 'coldim'};
    prod = prod.setup();
    
    divfv_T = SparseTensor();
    divfv_T = divfv_T.setFromTensorProd(normals, prod);
    divfv = divfv_T.getMatrix();
    
    %% We setup the consistent divergence operator
    % It consists of two parts,
    % divconsnf : nodefacecoltbl -> celltbl
    % divconsc : celltbl -> celltbl        
    
    % We agregate the contribution at each cell corner.
    % We use equal weights mcoef = (1/(number of nodes per cell)*(volume of the cell)).
    
    mcoef = 1./nnodespercell;    
    
    cno = celltbl.get('cells');
    vols = G.cells.volumes(cno);
    
    mcoef = vols/nnodespercell;
    
    prod = TensorProd();
    prod.tbl1 = celltbl;
    prod.tbl2 = cellnodefacecoltbl;
    prod.tbl3 = cellnodefacecoltbl;
    prod.mergefds = {'cells'};
    prod = prod.setup();
    
    mg = prod.eval(mcoef, g);
    
    % Multiply with Biot's coefficient alpha
    prod = TensorProd();
    prod.tbl1 = celltbl;
    prod.tbl2 = cellnodefacecoltbl;
    prod.tbl3 = cellnodefacecoltbl;
    prod.mergefds = {'cells'};
    prod = prod.setup();
    
    mg = prod.eval(alpha, mg);
    
    % We assemble divconsnf : nodefacecoltbl -> celltbl
    prod = TensorProd();
    prod.tbl1 = cellnodefacecoltbl; 
    prod.tbl2 = nodefacecoltbl;
    prod.tbl3 = celltbl;
    prod.reducefds = {'nodes', 'faces', 'coldim'};
    prod = prod.setup();
    
    divconsnf_T = SparseTensor();
    divconsnf_T = divconsnf_T.setFromTensorProd(mg, prod);
    divconsnf = divconsnf_T.getMatrix();
    
    % We assemble divconsc : cellcoltbl -> celltbl    
    prod = TensorProd();
    prod.tbl1 = cellnodefacecoltbl; 
    prod.tbl2 = cellcoltbl;
    prod.tbl3 = celltbl;
    prod.mergefds = {'cells'};
    prod.reducefds = {'coldim'};
    prod = prod.setup();
    
    divconsc_T = SparseTensor();
    % Beware of minus sign below
    divconsc_T = divconsc_T.setFromTensorProd(-mg, prod);
    divconsc = divconsc_T.getMatrix();
    
    assembly = struct('divfv'    , divfv    , ...
                      'divconsnf', divconsnf, ...
                      'divconsc' , divconsc);
    
end

