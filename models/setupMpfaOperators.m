function operators = setupMpfaOperators(model)
    
    G = model.G;
    mech = model.mech;
    prop = mech.prop;
    loadstruct = mech.loadstruct;
    
    eta = 0;
    bcetazero = false;
    
    [tbls, mappings] = setupStandardTables(G);

    assembly = assembleMPFA(G, K, bcstruct, src, eta, tbls, mappings, 'adoperators', true);
    
    operators = assembly.adoperators;
    
    tpfaoperators = setupOperatorsTPFA(G, rock);
    
    operators.pv = tpfaoperators.pv;
    
end

