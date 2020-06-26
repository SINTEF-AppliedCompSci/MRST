function operators = setupMpsaOperators(model)
    
    G = model.G;
    mech = model.mech;
    prop = mech.prop;
    loadstruct = mech.loadstruct;
    
    eta = 0;
    bcetazero = false;
    
    [tbls, mappings] = setupStandardTables(G);

    assembly = assembleMPSA(G, prop, loadstruct, eta, tbls, mappings, 'bcetazero', bcetazero, 'addAdOperators', true); 
    
    operators = assembly.adoperators;
    
end

