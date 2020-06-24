function operators = setupMpfaAdOperators(model)    
    
    %% We setup the usual operators and the mpfa flux operator
    %
    % we set up the mappings
    %
    %   F1 : nodefacetbl -> intfacetbl 
    %   F2 : celltbl -> intfacetbl
    %
    %  such that the flux u in intfacetbl (interior faces) is given by
    %
    %  u  = [F1  F2] * [ pnf (pressure at nodefacetbl);
    %                    pc  (pressure at celltbl)];
    %
    %  Then, we proceed with the reduction to remove dependency in pnf (pressure at nodefacetbl)
    %
    
    G         = model.G;
    rock      = model.rock;
    eta       = model.eta;
    bcetazero = model.bcetazero;
    
    perm = rock.perm;
    
    assert(numel(perm) == G.cells.num, 'only isotropic perm for the moment');
    
    [tbls, mappings] = setupStandardTables(G, 'useVirtual', false);

    celltbl = tbls.celltbl;
    colrowtbl = tbls.colrowtbl;
    cellcolrowtbl = tbls.cellcolrowtbl;%
    
    prod = TensorProd();
    prod.tbl1 = colrowtbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = cellcolrowtbl;
    prod = prod.setup();

    K = prod.eval([1; 0; 0; 1], perm);
    src = []; % no source at this stage
    bcstruct = []; % zero neumann bc
    
    assembly = assembleMPFA(G, K, bcstruct, src, eta, tbls, mappings, 'onlyAssemblyMatrices', true);
        
    matrices = assembly.matrices;
    nKg = assembly.nKg;
    
    invA11 = matrices.invA11;
    A12    = matrices.A12;
    
    cellfacetbl      = tbls.cellfacetbl;
    facetbl          = tbls.facetbl;
    cellnodeface2tbl = tbls.cellnodeface2tbl;
    nodefacetbl      = tbls.nodefacetbl;

    intfaces = find(all(G.faces.neighbors, 2));
    intfacetbl.faces = intfaces;
    intfacetbl = IndexArray(intfacetbl);
    
    map = TensorMap();
    map.fromTbl  = cellfacetbl;
    map.toTbl    = facetbl;
    map.mergefds = {'faces'};
    map = map.setup();
    
    ncellperface = map.eval(ones(cellfacetbl.num, 1));
    
    cfno = cellfacetbl.get('faces');
    sgn = 2*(cfno == G.faces.neighbors(cfno, 1)) - 1;

    prod = TensorProd();
    prod.tbl1 = facetbl;
    prod.tbl2 = cellfacetbl;
    prod.tbl3 = cellfacetbl;
    prod.mergefds = {'faces'};
    prod = prod.setup();
    
    wsgn = prod.eval(1./ncellperface, sgn);
    
    prod = TensorProd();
    prod.tbl1 = cellfacetbl;
    prod.tbl2 = cellnodeface2tbl;
    prod.tbl3 = cellnodeface2tbl;
    prod.replacefds1 = {{'faces', 'faces1'}};
    prod.mergefds = {'cells', 'faces1'};
    prod = prod.setup();
    
    wnKg = prod.eval(wsgn, nKg);

    gen = CrossIndexArrayGenerator();
    gen.tbl1 = cellnodeface2tbl;
    gen.tbl2 = intfacetbl;
    gen.replacefds2 = {{'faces', 'faces1'}};
    gen.mergefds = {'faces1'};
    
    cellnodeintface2tbl = gen.eval();
 
    map = TensorMap();
    map.fromTbl = cellnodeface2tbl;
    map.toTbl = cellnodeintface2tbl;
    map.mergefds = {'cells', 'nodes', 'faces1', 'faces2'};
    map = map.setup();
    
    wnKg = map.eval(wnKg);
    
    prod = TensorProd();
    prod.tbl1 = cellnodeintface2tbl;
    prod.tbl2 = nodefacetbl;
    prod.tbl3 = intfacetbl;
    prod.replacefds1 = {{'faces1', 'faces'}};
    prod.replacefds2 = {{'faces', 'faces2'}};
    prod.reducefds = {'nodes', 'faces2'};
    prod = prod.setup();
    
    F1_T = SparseTensor();
    F1_T = F1_T.setFromTensorProd(wnKg, prod);
    
    F1 = F1_T.getMatrix();

    prod = TensorProd();
    prod.tbl1 = cellnodeintface2tbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = intfacetbl;
    prod.replacefds1 = {{'faces1', 'faces'}};
    prod.reducefds = {'cells'};
    prod = prod.setup();
    
    F2_T = SparseTensor();
    % note the minus sign
    F2_T = F2_T.setFromTensorProd(-wnKg, prod);
    
    F2 = F2_T.getMatrix();
    
    %% We set up flux operator 
    %
    %  F : celltbl -> facetbl
    %
    % We assume Neumann boundary condition for flow so that we have
    %
    % [A11, A12] * [ pnf (pressure at nodefacetbl);
    %                pc  (pressure at celltbl)     ]   =   0;  
    %
    
    F = F2 - F1*invA11*A12;
    
    operators = setupOperatorsTPFA(G, rock);
    operators.mpfaKgrad = F;
    
end

