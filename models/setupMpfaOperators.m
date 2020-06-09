function operators = setupMpfaOperators(model)
    
    G = model.G;
    bcstruct = model.fluid.bcstruct;
    
    eta = 0;
    bcetazero = false;

    perm = model.rock.perm;
    
    assert(numel(perm) == G.cells.num, 'only isotropic perm for the moment');
    [tbls, mappings] = setupStandardTables(G, 'useVirtual', false);

    celltbl = tbls.celltbl;
    colrowtbl = tbls.colrowtbl;
    cellcolrowtbl = tbls.cellcolrowtbl;
    
    prod = TensorProd();
    prod.tbl1 = colrowtbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = cellcolrowtbl;
    prod = prod.setup();

    K = prod.eval([1; 0; 0; 1], perm);
    src = []; % no source at this stage
    
    assembly = assembleMPFA(G, K, bcstruct, src, eta, tbls, mappings, 'adoperators', true);
    
    operators = assembly.adoperators;
    
    % recover pore volume from TPFA (could be optimized!)
    rock = model.rock;
    tpfaoperators = setupOperatorsTPFA(G, rock);
    operators.pv = tpfaoperators.pv;
    
end

