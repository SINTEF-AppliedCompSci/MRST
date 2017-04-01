function [krW , krO , krG , krS , ...
          muW , muO , muG , muS , ...
          rhoW, rhoO, rhoG, rhoS, ...
          bW  , bO  , bG  , bS  , ...
          bW0 , bO0 , bG0 , bS0 , ...
          pvMult, transMult, mobMult, pvMult0, T] ...
                  =  getDynamicQuantitiesSolvent(model, p0, p, sW, sO, sG, sS)

    G     = model.G;
    fluid = model.fluid;
    op    = model.operators;

    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
    [krW, krO, krG, krS, sWres, sOres, sSGres] = computeRelPermSolvent(sW, sO, sG, sS, p, fluid, mobMult);
    
    % Compute transmissibility
    T = op.T.*transMult;

    [muW, muO, muG, muS, rhoW, rhoO, rhoG, rhoS] ...
        = computeViscositiesAndDensities(sW, sO, sG, sS, sWres, sOres, sSGres, fluid, p);
    
    bW = fluid.bW(p);
    bO = fluid.bO(p);
    bG = fluid.bG(p);
    bS = fluid.bS(p);
    
    bW0 = fluid.bW(p0);
    bO0 = fluid.bO(p0);
    bG0 = fluid.bG(p0);
    bS0 = fluid.bS(p0);
    
    
end

function [krW, krO, krG, krS, sWres, sOres, sSGres] = computeRelPermSolvent(sW, sO, sG, sS, p, fluid, mobMult)

    %% 
    
    krW = @(sW) fluid.krW(sW);
    
    krO_i = @(sO) fluid.krO(sO);
    krO_m = @(sO, sG, sS) sO./(sO + sG + sS).*fluid.krOW(sO + sG + sS);
    
    
    krGT_i = @(sG, sS) fluid.krG(sG + sS);
    krGT_m = @(sO, sG, sS) (sS + sG)./(sO + sG + sS).*fluid.krOW(sO + sG + sS);
    
    
    krG_i = @(sG, sS) krGT_i(sG, sS).*(sG./(sS + sG));
    krG_m = @(sO, sG, sS) krGT_m(sO, sG, sS).*(sG./(sS + sG));
    
    krS_i = @(sG, sS) krGT_i(sG, sS).*(sS./(sS + sG));
    krS_m = @(sO, sG, sS) krGT_m(sO, sG, sS).*(sS./(sS + sG));
    
    p = p.val;
    
    sW = sW.val; sO = sO.val; sS = sS.val;
    
    if isempty(sG)
        sG = zeros(numel(sW),1);
    else
        sG = sG.val;
    end
    
    M = fluid.Msat(sG, sS).*fluid.Mpres(p);

    sOres_m    = fluid.sOres_m ;    % Residual oil saturation     without solvent
    sOres_i    = fluid.sOres_i ;    % Residual oil saturation     without solvent
    sSGres_m   = fluid.sSGres_m;    % Residual oil saturation     without solvent
    sSGres_i   = fluid.sSGres_i;    % Residual oil saturation     without solvent
    sWres = fluid.sWres;
    
    % Interpolated water/oil residual saturations
    sOres  = M.*sOres_m + (1 - M).*sOres_i;
    sSGres = M.*sSGres_m + (1 - M).*sSGres_i;
    
    

    
    sW_eff = (sW - sWres )./(1 - sOres - sSGres - sWres );
    sO_eff = (sO - sOres )./(1 - sWres - sSGres - sOres );
    sG_eff = (sG - sSGres)./(1 - sWres - sOres  - sSGres);
    sS_eff = (sS - sSGres)./(1 - sWres - sOres  - sSGres);
    
    krW = krW(sW_eff);
    krO = M.*krO_m(sO_eff, sG_eff, sS_eff) + (1-M).*krO_i(sO_eff);
    krG = M.*krG_m(sO_eff, sG_eff, sS_eff) + (1-M).*krG_i(sG_eff, sS_eff);
    krS = M.*krS_m(sO_eff, sG_eff, sS_eff) + (1-M).*krS_i(sG_eff, sS_eff);
    
    % Modifiy relperm by mobility multiplier (if any)
    krW = mobMult.*krW;
    krO = mobMult.*krO;
    krG = mobMult.*krG;
    krS = mobMult.*krS;

end

function [muW_eff, muO_eff, muG_eff, muS_eff, rhoW_eff, rhoO_eff, rhoG_eff, rhoS_eff] ...
          = computeViscositiesAndDensities(sW, sO, sG, sS, sWres, sOres, sSGres, fluid, p)

    sW = sW.val;
    sO = sO.val;
    if isempty(sG)
        sG = zeros(numel(sO),1);
    else
        sG = sG.val;
    end
    sS = sS.val;
    
    
      
    muW = fluid.muW(p); muW = muW.val;
    muO = fluid.muO(p); muO = muO.val;
    muG = fluid.muG(p); muG = muG.val;
    muS = fluid.muS(p); muS = muS.val;
    
    sOn = sO - sOres;
    sGn = sG - sSGres;
    sSn = sS - sSGres;
    sNn = sOn + sGn + sSn;
    sOSn = sOn + sSn;
    sSGn = sSn + sGn;
    
    a = 1/4;
    muMOS = muO.*muS./(sSn./sOSn.*muS.^a + sSn./sOSn.*muO.^a).^4;
    muMSG = muS.*muG./(sSn./sSGn.*muG.^a + sGn./sSGn.*muS.^a).^4;
    muM   = muO.*muS.*muG./(sOn./sNn.*muS.^a.*muG.^a ...
                        + sSn./sNn.*muO.^a.*muG.^a ...
                        + sGn./sNn.*muO.^a.*muS.^a).^4;
                 
    omega = fluid.mixPar;
    muW_eff = muW;
    muO_eff = muO.^(1-omega).*muMOS.^omega;
    muG_eff = muG.^(1-omega).*muMSG.^omega;
    muS_eff = muS.^(1-omega).*muM.^omega;
    
    rhoW = fluid.bW(p).*fluid.rhoWS;
    rhoO = fluid.bO(p).*fluid.rhoOS;
    rhoG = fluid.bG(p).*fluid.rhoGS;
    rhoS = fluid.bS(p).*fluid.rhoSS;
    
    rhoW_eff = rhoW;
    
    same_visc = muO == muS | muG == muS;

    [rhoO_eff, rhoG_eff, rhoS_eff] = deal(zeros(numel(sO),1));
    
    sN = sO + sG + sS;
    
    rhoM = rhoO.*sO./sN + rhoG.*sG./sN + rhoS.*sS./sN;
    rhoO_eff(same_visc) = (1-omega)*rhoO(same_visc) + omega*rhoM(same_visc);
    rhoG_eff(same_visc) = (1-omega)*rhoG(same_visc) + omega*rhoM(same_visc);
    rhoS_eff(same_visc) = (1-omega)*rhoS(same_visc) + omega*rhoM(same_visc);

    a = 1/4;
    sOsN_oe = muO.^a.*(muO_eff.^a - muS.^a)./(muO_eff.^a.*(muO.^a - muS.^a));
    sOsN_ge = muS.^a.*(muG_eff.^a - muG.^a)./(muG_eff.^a.*(muS.^a - muG.^a));
    sOGn = sOn + sGn;
    sOf = sOn./sOGn;
    sGf = sGn./sOGn;
    sSsN_se = (muS.^a.*(sGf.*muO.^a + sOf.*muG.^a) - muO.^a.*muG.^a.*(muS./muS_eff).^a)...
            ./(muS.^a.*(sGf.*muO.^a + sOf.*muG.^a) - muO.^a.*muG.^a);

    rhoO_eff = sOsN_oe(~same_visc).*rhoO(~same_visc) + (1-sOsN_oe(~same_visc)).*rhoS(~same_visc);
    rhoG_eff = sOsN_ge(~same_visc).*rhoS(~same_visc) + (1-sOsN_ge(~same_visc)).*rhoG(~same_visc);
    rhoS_eff = sSsN_se(~same_visc).*rhoS(~same_visc) + (1-sSsN_se(~same_visc)).*(rhoG(~same_visc).*sGf(~same_visc) + rhoO(~same_visc).*sOf(~same_visc));


end
    
    
    
    
    
    