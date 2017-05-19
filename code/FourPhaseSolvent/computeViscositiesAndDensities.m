function [muW_eff, muO_eff, muG_eff, muS_eff, rhoW_eff, rhoO_eff, rhoG_eff, rhoS_eff] = computeViscositiesAndDensities(fluid, p , sO , sG , sS , sOres , sSGres )
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
    
%     % Add small value to avoid 0/0-type expressions
%     tol = 10*eps;
%     sOn = sOn + (abs(sOn) < tol).*tol;
%     sGn = sGn + (abs(sGn) < tol).*tol;
%     sSn = sSn + (abs(sSn) < tol).*tol;
    
    sNn = sOn + sGn + sSn;
    sOSn = sOn + sSn;
    sSGn = sSn + sGn;
    
    % Calculate mixed viscosities
%     a = 1/4;
%     muMOS = muO.*muS./(sOn./sOSn.*muS.^a + sSn./sOSn.*muO.^a).^4;
%     muMSG = muS.*muG./(sSn./sSGn.*muG.^a + sGn./sSGn.*muS.^a).^4;
%     muM   = muO.*muS.*muG./(sOn./sNn.*muS.^a.*muG.^a ...
%                           + sSn./sNn.*muO.^a.*muG.^a ...
%                           + sGn./sNn.*muO.^a.*muS.^a).^4;
%       
    sOnsOSn = saturationFraction(sOn, sOSn);
    sSnsSGn = saturationFraction(sSn, sSGn);
    sOnsNn = saturationFraction(sOn, sNn);
    sSnsNn = saturationFraction(sSn, sNn);
    
    % Calculate mixed viscosities
    a = 1/4;
    muMOS = muO.*muS./(sOnsOSn.*muS.^a + (1 - sOnsOSn).*muO.^a).^4;
    muMSG = muS.*muG./(sSnsSGn.*muG.^a + (1 - sSnsSGn).*muS.^a).^4;
    muM   = muO.*muS.*muG./(sOnsNn.*muS.^a.*muG.^a ...
                          + sSnsNn.*muO.^a.*muG.^a ...
                          + (1 - (sOnsNn + sSnsNn)).*muO.^a.*muS.^a).^4;
  
    
    
    omega = fluid.mixPar;
    muW_eff = muW;
    muO_eff = muO.^(1-omega).*muMOS.^omega;
    muG_eff = muG.^(1-omega).*muMSG.^omega;
    muS_eff = muS.^(1-omega).*muM.^omega;
    
    % Unmixed densities at reservoir conditions
    rhoW = fluid.bW(p).*fluid.rhoWS;
    rhoO = fluid.bO(p).*fluid.rhoOS;
    rhoG = fluid.bG(p).*fluid.rhoGS;
    rhoS = fluid.bS(p).*fluid.rhoSS;
    
    % Effective fractional saturations
    sOsN_Oeff = muO.^a.*(muO_eff.^a - muS.^a)./(muO_eff.^a.*(muO.^a - muS.^a));
    sOsN_Geff = muS.^a.*(muG_eff.^a - muG.^a)./(muG_eff.^a.*(muS.^a - muG.^a));
%     sOGn      = sOn + sGn;
%     sOf       = sOn./sOGn;
%     sGf       = sGn./sOGn;
%     sSsN_Seff = (muS.^a.*(sGf.*muO.^a + sOf.*muG.^a) - muO.^a.*muG.^a.*(muS./muS_eff).^a)...
%                   ./(muS.^a.*(sGf.*muO.^a + sOf.*muG.^a) - muO.^a.*muG.^a);

    sOGn = sOn + sGn;
    tol = eps;
    if all(sGn < tol)
        sGf = 0;
    else
        sGf = saturationFraction(sGn, sOGn);
    end
    
    sSsN_Seff = (muS.^a.*(sGf.*muO.^a + (1-sGf).*muG.^a) - muO.^a.*muG.^a.*(muS./muS_eff).^a)...
                  ./(muS.^a.*(sGf.*muO.^a + (1-sGf).*muG.^a) - muO.^a.*muG.^a);

    % Expressions are sinuglar if muO == muG, in which case we replace the
    % by a simple interpolation rho*(1-omega) + rhoM*omega
    tol = 1e-10;
    eq = abs(muO - muS) < tol | abs(muS - muG) < tol;
    
%     sO = sO + (abs(sO) < tol).*tol;
%     sG = sG + (abs(sG) < tol).*tol;
%     sS = sS + (abs(sS) < tol).*tol;
%     rhoM = rhoO.*sO./sN + rhoG.*sG./sN + rhoS.*sS./sN;
    sN = sO + sG + sS;

    sOsN = saturationFraction(sO, sN);
    sGsN = saturationFraction(sG, sN);
    rhoM = rhoO.*sOsN + rhoG.*sGsN + rhoS.*(1 - (sOsN + sGsN));
    
    % Calulcate mixed densities
    rhoW_eff = rhoW;
    rhoO_eff = (sOsN_Oeff.*rhoO + (1-sOsN_Oeff).*rhoS).*(~eq) ...
                                       + ((1-omega)*rhoO + omega*rhoM).*eq;
    rhoG_eff = (sOsN_Geff.*rhoS + (1-sOsN_Geff).*rhoG).*(~eq) ...
                                       + ((1-omega)*rhoG + omega*rhoM).*eq;
    rhoS_eff = (sSsN_Seff.*rhoS + (1-sSsN_Seff).*(rhoG.*sGf + rhoO.*(1-sGf))).*(~eq) ...
                                       + ((1-omega)*rhoS + omega*rhoM).*eq;
    
end

function f = saturationFraction(sNum, sDen)

    tol = eps;
    f = sNum./(sDen + (sDen < tol).*tol);
    
%     num = sNum - sDen;
%     f = 1 - (num + (num < tol).*tol)./(sDen + (sDen < tol).*tol);
    
    
end