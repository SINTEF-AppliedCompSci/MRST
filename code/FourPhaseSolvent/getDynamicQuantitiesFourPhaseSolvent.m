function [kr_eff, mu_eff, rho_eff, b_eff, b0_eff, pvMult, pvMult0, T] ...
    =  getDynamicQuantitiesFourPhaseSolvent(model, p, p0, sW, sO, sG, sS, sW0, sO0, sG0, sS0)

    fluid = model.fluid;
    op    = model.operators;

    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
    T = op.T.*transMult;
    
    [sWres, sOres, sSGres] = residualSaturations(fluid, p, sG, sS);
    [~    , sO0res, sSG0res] = residualSaturations(fluid, p0, sG0, sS0);
    
    kr_eff = computeRelPermSolvent(fluid, p, sW, sO, sG, sS, sWres, sOres, sSGres, mobMult);
    
    [mu_eff, rho_eff] ...
        = computeViscositiesAndDensities(fluid, p, sO, sG, sS, sOres, sSGres);
    
    [~     , rho0_eff] ...
        = computeViscositiesAndDensities(fluid, p0, sO0, sG0, sS0, sO0res, sSG0res);
    
    b_eff = calculateFormationVolumeFactors(fluid, p, rho_eff);
    b0_eff = calculateFormationVolumeFactors(fluid, p0, rho0_eff);
    
end

function [sWres, sOres, sSGres] = residualSaturations(fluid, p, sG, sS)

     % Residual saturations for the immiscible and miscible extrema
    sOres_m    = fluid.sOres_m ;
    sOres_i    = fluid.sOres_i ;
    sSGres_m   = fluid.sSGres_m;
    sSGres_i   = fluid.sSGres_i;
    
    M = fluid.Msat(sG, sS).*fluid.Mpres(p);
    
    if isprop(M, 'val')
        M.val(isnan(M.val)) = 1;
    else
        M(isnan(M)) = 1;
    end
    
    % Interpolated water/oil residual saturations
    sWres = fluid.sWres;
    sOres  = M.*sOres_m + (1 - M).*sOres_i;
    sSGres = M.*sSGres_m + (1 - M).*sSGres_i;
    
end

function kr_eff = computeRelPermSolvent(fluid, p, sW, sO, sG, sS, sWres, sOres, sSGres, mobMult)

    %% 
    
    % Water relperm not affected by the solvent
    krW = @(sW) fluid.krW(sW);
    
    % Relperm of hydrocarbon to water as afunciton of total HC saturation
    krN = @(sO, sG, sS) fluid.krOW(sO + sG + sS);
    
    % Relperm of oil, immiscible (i) and miscible (m)
    krO_i = @(sO) fluid.krO(sO);
    krO_m = @(sO, sG, sS) sO./(sO + sG + sS).*krN(sO, sG, sS);
    
    % Total relperm og gas and solvent
    krGT_i = @(sG, sS) fluid.krG(sG + sS);
    krGT_m = @(sO, sG, sS) (sG + sS)./(sO + sG + sS).*krN(sO, sG, sS);
    
    % Relperm of gas
    krG_i = @(sG, sS) krGT_i(sG, sS).*(sG./(sS + sG));
    krG_m = @(sO, sG, sS) krGT_m(sO, sG, sS).*(sG./(sS + sG));
    
    % Relperm of solvent
    krS_i = @(sG, sS) krGT_i(sG, sS).*(sS./(sS + sG));
    krS_m = @(sO, sG, sS) krGT_m(sO, sG, sS).*(sS./(sS + sG));
        
    sW_eff = max((sW - sWres )./(1 - sOres - sSGres - sWres ),0);
    sO_eff = max((sO - sOres )./(1 - sWres - sSGres - sOres ),0);
    sG_eff = max((sG - sSGres)./(1 - sWres - sOres  - sSGres),0);
    sS_eff = max((sS - sSGres)./(1 - sWres - sOres  - sSGres),0);
    
    M = fluid.Msat(sG, sS).*fluid.Mpres(p);
    M.val(isnan(M.val)) = 1;

    
    krW_eff = krW(sW_eff);
    krO_eff = M.*krO_m(sO_eff, sG_eff, sS_eff) + (1-M).*krO_i(sO_eff);
    krG_eff = M.*krG_m(sO_eff, sG_eff, sS_eff) + (1-M).*krG_i(sG_eff, sS_eff);
    krS_eff = M.*krS_m(sO_eff, sG_eff, sS_eff) + (1-M).*krS_i(sG_eff, sS_eff);
    
    krW_eff.val(isnan(krW_eff.val)) = 1;
    krO_eff.val(isnan(krO_eff.val)) = 1;
    krG_eff.val(isnan(krG_eff.val)) = 1;
    krS_eff.val(isnan(krS_eff.val)) = 1;
    
    % Modifiy relperm by mobility multiplier (if any)
    krW_eff = mobMult.*krW_eff;
    krO_eff = mobMult.*krO_eff;
    krG_eff = mobMult.*krG_eff;
    krS_eff = mobMult.*krS_eff;
    
    kr_eff = {krW_eff, krO_eff, krG_eff, krS_eff};

end

function [mu_eff, rho_eff] ...
    = computeViscositiesAndDensities(fluid, p, sO, sG, sS, sOres, sSGres)
    
    muW = fluid.muW(p);
    muO = fluid.muO(p);
    muG = fluid.muG(p);
    muS = fluid.muS(p);
    
    sOn = max(sO - sOres,0);
    sGn = max(sG - sSGres,0);
    sSn = max(sS - sSGres,0);
    sNn = sOn + sGn + sSn;
    sOSn = sOn + sSn;
    sSGn = sSn + sGn;
    
    a = 1/4;
    muMOS = muO.*muS./(sOn./sOSn.*muS.^a + sSn./sOSn.*muO.^a).^4;
    muMSG = muS.*muG./(sSn./sSGn.*muG.^a + sGn./sSGn.*muS.^a).^4;
    muM   = muO.*muS.*muG./(sOn./sNn.*muS.^a.*muG.^a ...
                        + sSn./sNn.*muO.^a.*muG.^a ...
                        + sGn./sNn.*muO.^a.*muS.^a).^4;
    
    if isprop(muMOS, 'val')
        muMOS.val(isnan(muMOS.val)) = 1;
        muMSG.val(isnan(muMSG.val)) = 1;
        muM.val(isnan(muM.val)) = 1;
    else
        muMOS(isnan(muMOS)) = 1;
        muMSG(isnan(muMSG)) = 1;
        muM(isnan(muM)) = 1;
    end
                 
    omega = fluid.mixPar;
    muW_eff = muW;
    muO_eff = muO.^(1-omega).*muMOS.^omega;
    muG_eff = muG.^(1-omega).*muMSG.^omega;
    muS_eff = muS.^(1-omega).*muM.^omega;
    
    mu_eff = {muW_eff, muO_eff, muG_eff, muS_eff};
    
    rhoW = fluid.bW(p).*fluid.rhoWS;
    rhoO = fluid.bO(p).*fluid.rhoOS;
    rhoG = fluid.bG(p).*fluid.rhoGS;
    rhoS = fluid.bS(p).*fluid.rhoSS;
    
    sN = sO + sG + sS;
    
    rhoM = rhoO.*sO./sN + rhoG.*sG./sN + rhoS.*sS./sN;

    a = 1/4;
    sOsN_Oeff = muO.^a.*(muO_eff.^a - muS.^a)./(muO_eff.^a.*(muO.^a - muS.^a));
    sOsN_Geff = muS.^a.*(muG_eff.^a - muG.^a)./(muG_eff.^a.*(muS.^a - muG.^a));
    sOGn = sOn + sGn;
    sOf = sOn./sOGn;
    sGf = sGn./sOGn;
    sSsN_se = (muS.^a.*(sGf.*muO.^a + sOf.*muG.^a) - muO.^a.*muG.^a.*(muS./muS_eff).^a)...
                  ./(muS.^a.*(sGf.*muO.^a + sOf.*muG.^a) - muO.^a.*muG.^a);

    % Expressions are sinuglar if muO == muG
    tol = 1e-10;
    eq = abs(muO - muS) < tol | abs(muS - muG) < tol;
    
    rhoW_eff = rhoW;
    rhoO_eff = (sOsN_Oeff.*rhoO + (1-sOsN_Oeff).*rhoS).*(~eq) ...
                                       + ((1-omega)*rhoO + omega*rhoM).*eq;
    rhoG_eff = (sOsN_Geff.*rhoS + (1-sOsN_Geff).*rhoG).*(~eq) ...
                                       + ((1-omega)*rhoG + omega*rhoM).*eq;
    rhoS_eff = (sSsN_se.*rhoS + (1-sSsN_se).*(rhoG.*sGf + rhoO.*sOf)).*(~eq) ...
                                       + ((1-omega)*rhoS + omega*rhoM).*eq;
                                   
    rho_eff = {rhoW_eff, rhoO_eff, rhoG_eff, rhoS_eff};
    
end
    
function b_eff = calculateFormationVolumeFactors(fluid, p, rho_eff)

    rhoO_eff = rho_eff{2};
    rhoG_eff = rho_eff{3};
    rhoS_eff = rho_eff{4};
    

    bW_eff = fluid.bW(p);
    bO_eff = rhoO_eff./fluid.rhoOS;
    bG_eff = rhoG_eff./fluid.rhoGS;
    bS_eff = rhoS_eff./fluid.rhoSS;
    
    
    b_eff = {bW_eff, bO_eff, bG_eff, bS_eff};
    
end
    