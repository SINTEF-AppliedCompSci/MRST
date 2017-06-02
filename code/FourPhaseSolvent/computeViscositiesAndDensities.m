function [muW_eff, muO_eff, muG_eff, muS_eff, rhoW_eff, rhoO_eff, rhoG_eff, rhoS_eff] = computeViscositiesAndDensities(fluid, p , sO , sG , sS , sOres , sSGres )
    % Calculates effective viscosities and densities using Todd-Longstaff
    % model + 1/4th-power mixing rule

    % Unmixed viscosities at reservoir conditions
    muW = fluid.muW(p);
    muO = fluid.muO(p);
    muG = fluid.muG(p);
    muS = fluid.muS(p);
    
    % Subtract residual sturations
    sOn = max(sO - sOres, 0);
    sGn = max(sG - sSGres, 0);
    sSn = max(sS - sSGres, 0);

    sSnsOSn = saturationFraction(sSn, sOn);
    sSnsSGn = saturationFraction(sSn, sGn);
    sSnsNn = saturationFraction(sSn, sOn + sGn);
    sGnsNn = saturationFraction(sGn, sOn + sSn);
    
    
    % Calculate mixed viscosities
    a = 1/4;
    muSa = muS.^a;
    muGa = muG.^a;
    muOa = muO.^a;
    muMOS = muO.*muS./((1-sSnsOSn).*muSa + sSnsOSn.*muOa).^4;
    muMSG = muS.*muG./(sSnsSGn.*muGa + (1-sSnsSGn).*muSa).^4;
    muM   = muO.*muS.*muG./((1-(sGnsNn + sSnsNn)).*muSa.*muGa ...
                          + sSnsNn.*muOa.*muGa ...
                          + sGnsNn.*muOa.*muSa).^4;
  
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
    
    muOmuS = (muO./muS).^a;
    muSmuO = 1./muOmuS;
    muSmuG = (muS./muG).^a;
    sOsN_Oeff = max((muOmuS - (muO./muO_eff).^a)./(muOmuS-1),0);
    sOsN_Geff = max((muSmuG - (muS./muG_eff).^a)./(muSmuG-1),0);
    
    sGf = saturationFraction(sGn, sOn);
    
    sSsN_Seff = max((muSmuG.*sGf + muSmuO.*(1-sGf) - (muS./muS_eff).^a)...
               ./(muSmuG.*sGf + muSmuO.*(1-sGf) - 1),0);
    
    % Expressions are sinuglar if muO == muG, in which case we replace the
    % by a simple interpolation rho*(1-omega) + rhoM*omega
    tol = 1e-10;
    eq = abs(muO - muS) < tol | abs(muS - muG) < tol;

    sOsN = saturationFraction(sO, sG + sS);
    sGsN = saturationFraction(sG, sO + sS);
    
    
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
% 
% function s = trimSaturations(s)
% 
%     tol = eps;
%     s = max(s + tol*(s<tol),0);
%     
% end

% function f = saturationFraction(sNum, sDen)
% 
% %     tol = 10*eps;
% %     zz = sNum < tol & sDen < tol;
% %     f = sNum./(sDen + (sDen < tol).*tol).*(~zz) + 0.*zz;
% 
%     tol = 10*eps;
%     zz = sNum < tol & sDen < tol;
%     sNum = sNum - (sNum < tol).*sNum;
%     sDen = sDen - (sDen < tol).*sDen;
%     f = sNum./(sDen + (sDen < tol).*tol).*(~zz) + 0.*zz;
%        
% end
