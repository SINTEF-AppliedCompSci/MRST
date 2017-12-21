function state0 = convertBlackOilStateToCompositional(bomodel, state)
% Convert BO state to compositional-like state
    [p, sO, sG, sW] = bomodel.getProps(state, 'pressure', 'sO', 'sG', 'sW');
    if bomodel.disgas
        rs = bomodel.getProp(state, 'rs');
    else
        rs = 0;
    end
    if bomodel.vapoil
        rv = bomodel.getProp(state, 'rv');
    else
        rv = 0;
    end
    [xo, xg, yo, yg, zo, zg, rhoO, rhoG] = blackOilToMassFraction(bomodel, p, sO, sG, rs, rv);
    
    state0 = initResSol(bomodel.G, p, [sW, sO, sG]);
    state0.components = [zo, zg];
    state0.x = [xo, xg];
    state0.L = rhoO.*sO./(rhoO.*sO + rhoG.*sG);
    
    state0.y = repmat([yo, yg], bomodel.G.cells.num, 1);
    state0.T = repmat(273.15, bomodel.G.cells.num, 1);
end