function [kr_eff, mu_eff, rho_eff, b_eff, b0_eff, pvMult, pvMult0, T] ...
    =  getDynamicQuantitiesFourPhaseSolvent(model, p, p0, sW, sO, sG, sS, sW0, sO0, sG0, sS0)
    % Calculates dynamic quantities (kr, mu, rho, b, transmissiblity and
    % multipliers.

    fluid = model.fluid;
    op    = model.operators;

    % Get multipliers
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
    T = op.T.*transMult;
    
    % Calculate residual saturations
    [sWres, sOres , sSGres ] = residualSaturations(fluid, p, sG, sS);
    [~    , sO0res, sSG0res] = residualSaturations(fluid, p0, sG0, sS0);
    
    % Calculate effective relperms
    kr_eff = computeRelPermSolvent(fluid, p, sW, sO, sG, sS, sWres, sOres, sSGres, mobMult);
    
    % Calulate effective viscosities and densities
    [mu_eff, rho_eff ] = computeViscositiesAndDensities(fluid, p , sO , sG , sS , sOres , sSGres );
    [~     , rho0_eff] = computeViscositiesAndDensities(fluid, p0, sO0, sG0, sS0, sO0res, sSG0res);
    
    % Calulcate effective formation volume factors
    b_eff  = computeFormationVolumeFactors(fluid, p , rho_eff );
    b0_eff = computeFormationVolumeFactors(fluid, p0, rho0_eff);
    
end

%% Helper functions

%--------------------------------------------------------------------------

function [sWres, sOres, sSGres] = residualSaturations(fluid, p, sG, sS)
    % Calculate effective residual saturations

    % Residual saturations for the immiscible and miscible extrema
    sOres_m    = fluid.sOres_m ;
    sOres_i    = fluid.sOres_i ;
    sSGres_m   = fluid.sSGres_m;
    sSGres_i   = fluid.sSGres_i;

    % Add small value to avoid 0/0-type expressions
    tol = 10*eps;
    sS = sS + (abs(sS) < tol).*tol;
    sG = sG + (abs(sG) < tol).*tol;
    
    % Misscibility is a funciton of the solvent fraction in the total gas
    % phase

    M = fluid.Msat(sG, sS).*fluid.Mpres(p);
    
    % Interpolated water/oil residual saturations
    sWres  = fluid.sWres;
    sOres  = M.*sOres_m  + (1 - M).*sOres_i ;
    sSGres = M.*sSGres_m + (1 - M).*sSGres_i;
    
end

%--------------------------------------------------------------------------

function kr_eff = computeRelPermSolvent(fluid, p, sW, sO, sG, sS, sWres, sOres, sSGres, mobMult)
    % Calulates effective relative permeabilities.
    
    % Relperm of hydrocarbon to water as a funciton of total HC saturation
    krN = fluid.krO(sO + sG + sS);
    
    % Immiscible relperms for oil and total gas
    krO_i = fluid.krO(sO);
    krGT_i = fluid.krG(sG + sS);
    
    % Add small value to avoid 0/0-type expressions
    tol = 10*eps;
    sO = sO + (abs(sO) < tol).*tol;
    sG = sG + (abs(sG) < tol).*tol;
    sS = sS + (abs(sS) < tol).*tol;
    
    % Immiscible gas and solvent relperms
    krG_i = sG./(sS + sG).*krGT_i;
    krS_i = sS./(sS + sG).*krGT_i;
    
    % Miscible relperms
    krO_m = sO./(sO + sG + sS).*krN;
    krGT_m = (sG + sS)./(sO + sG + sS).*krN;
    krG_m = sG./(sS + sG).*krGT_m;
    krS_m = sS./(sS + sG).*krGT_m;
    
    % Interpolate between miscible and immiscible cases (water relperm
    % not affected by the solvent)
    M = fluid.Msat(sG, sS).*fluid.Mpres(p);
    krW_eff = fluid.krW(sW);
    krO_eff = M.*krO_m + (1-M).*krO_i;
    krG_eff = M.*krG_m + (1-M).*krG_i;
    krS_eff = M.*krS_m + (1-M).*krS_i;
    
    % Modifiy relperm by mobility multiplier (if any)
    krW_eff = mobMult.*krW_eff;
    krO_eff = mobMult.*krO_eff;
    krG_eff = mobMult.*krG_eff;
    krS_eff = mobMult.*krS_eff;
    
    kr_eff = {krW_eff, krO_eff, krG_eff, krS_eff};

end

%--------------------------------------------------------------------------

function [mu_eff, rho_eff] ...
    = computeViscositiesAndDensities(fluid, p, sO, sG, sS, sOres, sSGres)
    % Calculates effective viscosities and densities using Todd-Longstaff
    % model + 1/4th-power mixing rule

    % Unmixed viscosities at reservoir conditions
    muW = fluid.muW(p);
    muO = fluid.muO(p);
    muG = fluid.muG(p);
    muS = fluid.muS(p);
    
    % Subtract residual sturations
    sOn = max(sO - sOres ,0);
    sGn = max(sG - sSGres,0);
    sSn = max(sS - sSGres,0);
    
    % Add small value to avoid 0/0-type expressions
    tol = 10*eps;
    sOn = sOn + (abs(sOn) < tol).*tol;
    sGn = sGn + (abs(sOn) < tol).*tol;
    sSn = sSn + (abs(sSn) < tol).*tol;
    
    sNn = sOn + sGn + sSn;
    sOSn = sOn + sSn;
    sSGn = sSn + sGn;
    
    % Calculate mixed viscosities
    a = 1/4;
    muMOS = muO.*muS./(sOn./sOSn.*muS.^a + sSn./sOSn.*muO.^a).^4;
    muMSG = muS.*muG./(sSn./sSGn.*muG.^a + sGn./sSGn.*muS.^a).^4;
    muM   = muO.*muS.*muG./(sOn./sNn.*muS.^a.*muG.^a ...
                          + sSn./sNn.*muO.^a.*muG.^a ...
                          + sGn./sNn.*muO.^a.*muS.^a).^4;
      
    omega = fluid.mixPar;
    muW_eff = muW;
    muO_eff = muO.^(1-omega).*muMOS.^omega;
    muG_eff = muG.^(1-omega).*muMSG.^omega;
    muS_eff = muS.^(1-omega).*muM.^omega;
    
    mu_eff = {muW_eff, muO_eff, muG_eff, muS_eff};
    
    % Unmixed densities at reservoir conditions
    rhoW = fluid.bW(p).*fluid.rhoWS;
    rhoO = fluid.bO(p).*fluid.rhoOS;
    rhoG = fluid.bG(p).*fluid.rhoGS;
    rhoS = fluid.bS(p).*fluid.rhoSS;
    
    % Effective fractional saturations
    sOsN_Oeff = muO.^a.*(muO_eff.^a - muS.^a)./(muO_eff.^a.*(muO.^a - muS.^a));
    sOsN_Geff = muS.^a.*(muG_eff.^a - muG.^a)./(muG_eff.^a.*(muS.^a - muG.^a));
    sOGn      = sOn + sGn;
    sOf       = sOn./sOGn;
    sGf       = sGn./sOGn;
    sSsN_Seff = (muS.^a.*(sGf.*muO.^a + sOf.*muG.^a) - muO.^a.*muG.^a.*(muS./muS_eff).^a)...
                  ./(muS.^a.*(sGf.*muO.^a + sOf.*muG.^a) - muO.^a.*muG.^a);

    % Expressions are sinuglar if muO == muG, in which case we replace the
    % by a simple interpolation rho*(1-omega) + rhoM*omega
    tol = 1e-10;
    eq = abs(muO - muS) < tol | abs(muS - muG) < tol;
    sN = sO + sG + sS;
    rhoM = rhoO.*sO./sN + rhoG.*sG./sN + rhoS.*sS./sN;
    
    % Calulcate mixed densities
    rhoW_eff = rhoW;
    rhoO_eff = (sOsN_Oeff.*rhoO + (1-sOsN_Oeff).*rhoS).*(~eq) ...
                                       + ((1-omega)*rhoO + omega*rhoM).*eq;
    rhoG_eff = (sOsN_Geff.*rhoS + (1-sOsN_Geff).*rhoG).*(~eq) ...
                                       + ((1-omega)*rhoG + omega*rhoM).*eq;
    rhoS_eff = (sSsN_Seff.*rhoS + (1-sSsN_Seff).*(rhoG.*sGf + rhoO.*sOf)).*(~eq) ...
                                       + ((1-omega)*rhoS + omega*rhoM).*eq;
                                   
    rho_eff = {rhoW_eff, rhoO_eff, rhoG_eff, rhoS_eff};
    
end

%--------------------------------------------------------------------------
    
function b_eff = computeFormationVolumeFactors(fluid, p, rho_eff)
    % Calculate effective (inverse) formation volume factors due to new densities

    % Effective densities
    rhoO_eff = rho_eff{2};
    rhoG_eff = rho_eff{3};
    rhoS_eff = rho_eff{4};
    
    % New formation volume factors b = rho_eff/rhoS
    bW_eff = fluid.bW(p);
    bO_eff = rhoO_eff./fluid.rhoOS;
    bG_eff = rhoG_eff./fluid.rhoGS;
    bS_eff = rhoS_eff./fluid.rhoSS;
    
    b_eff = {bW_eff, bO_eff, bG_eff, bS_eff};
    
end
    