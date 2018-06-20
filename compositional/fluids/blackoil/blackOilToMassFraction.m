function [xo, xg, yo, yg, zo, zg, rhoO, rhoG, muO, muG, freeGas] = blackOilToMassFraction(model, p, sO, sG, rs, rv)
% Evaluate properties
    fluid = model.fluid;
    rhoGS = fluid.rhoGS;
    rhoOS = fluid.rhoOS;
    
    disgas = model.disgas;
    if disgas
        isSatRs = rs >= fluid.rsSat(p);
        bO = fluid.bO(p, rs, isSatRs);

        muO = fluid.muO(p, rs, isSatRs);
        rhoO = bO.*(rs*rhoGS + rhoOS);
    else
        bO = fluid.bO(p);
        muO = fluid.muO(p);
        rhoO = bO.*rhoOS;
    end

    bG = fluid.bG(p);
    muG = fluid.muG(p);
    
    rhoG = bG.*rhoGS;
    zg = (rhoG.*sG + disgas.*rhoGS.*rs.*bO.*sO)./(rhoO.*sO + rhoG.*sG);
    zo = 1 - zg;
    
    xo = bO.*rhoOS./rhoO;
    xg = disgas.*bO.*rhoGS.*rs./rhoO;
    yg = 1;
    yo = 0;
    
    if model.disgas
        freeGas = isSatRs;
    else
        freeGas = zg > 0;
    end
end