function [muW_eff, muO_eff, muG_eff, muS_eff, rhoW_eff, rhoO_eff, rhoG_eff, rhoS_eff] = computeViscositiesAndDensities(fluid, p , sO , sG , sS , sOres , sSGres )
    % Calculates effective viscosities and densities using Todd-Longstaff
    % model + 1/4th-power mixing rule

    % Unmixed viscosities at reservoir conditions
    muW = fluid.muW(p);
    muO = fluid.muO(p);
    muG = fluid.muG(p);
    muS = fluid.muS(p);
    
    tol = 1e-10;
    sO(sO < tol) = 0;
    sG(sG < tol) = 0;
    sS(sS < tol) = 0; 
    
    % Caluculate mobile saturations
    sOn = max(sO - sOres, 0);
    sGn = max(sG - sSGres, 0);
    sSn = max(sS - sSGres, 0);
    
    
    % Calculate mixed viscosities
    a = 1/4;
    muSa = muS.^a;
    muGa = muG.^a;
    muOa = muO.^a;

    sSnsOSn = sSn./(sOn + sSn);
    sSnsOSn(isnan(double(sSnsOSn))) = 0;
    sSnsSGn = sSn./(sGn + sSn);
    sSnsSGn(isnan(double(sSnsSGn))) = 0;
    sOnsNn = sOn./(sOn + sGn + sSn);
    sOnsNn(isnan(double(sOnsNn))) = 0;
    sGnsNn = sGn./(sOn + sGn + sSn);
    sGnsNn(isnan(double(sGnsNn))) = 0;
        
    muMOS = muO.*muS./((1-sSnsOSn).*muSa + sSnsOSn.*muOa).^4;
    muMSG = muS.*muG./(sSnsSGn.*muGa + (1-sSnsSGn).*muSa).^4;
    muM   = muO.*muS.*muG./(sOnsNn.*muSa.*muGa ...
                          + (1-sOnsNn-sGnsNn).*muOa.*muGa ...
                          + sGnsNn.*muOa.*muSa).^4;
    
    % Effective viscosities are determined by the mixing parameter
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

    sGf = sGn./(sOn + sGn);
    sGf(isnan(double(sGf))) = 0;

    sSsN_Seff = max((muSmuG.*sGf + muSmuO.*(1-sGf) - (muS./muS_eff).^a)...
               ./(muSmuG.*sGf + muSmuO.*(1-sGf) - 1),0);
           
    % Expressions are sinuglar if muO == muG, in which case we replace the
    % by a simple interpolation rho*(1-omega) + rhoM*omega
    eq = abs(muO - muS) < 1e-10 | abs(muS - muG) < 1e-10;

    sOsN = sO./(sO + sG + sS);
    sOsN(isnan(double(sOsN))) = 0;
    sGsN = sG./(sO + sG + sS);
    sGsN(isnan(double(sGsN))) = 0;
    
    rhoM = rhoO.*sOsN + rhoG.*sGsN + rhoS.*(1 - (sOsN + sGsN));
    
    % Calulcate mixed densities
    rhoW_eff = rhoW;
    
    rhoO_eff = (sOsN_Oeff.*rhoO + (1-sOsN_Oeff).*rhoS).*(~eq) ...
                                       + ((1-omega)*rhoO + omega*rhoM).*eq;
                                   
    rhoG_eff = (sOsN_Geff.*rhoS + (1-sOsN_Geff).*rhoG).*(~eq) ...
                                       + ((1-omega)*rhoG + omega*rhoM).*eq;
    
    sGf = sGn./(sOn + sGn);
    sGf(isnan(double(sGf))) = 0;
    
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

%     tol = 10*eps;
%     sOn(sOn < tol) = 0;
%     sGn(sGn < tol) = 0;
%     sSn(sSn < tol) = 0;
    
%     tol = 1e-2;
% % %     tol = eps;
%     sSnsOSn = saturationFraction(sSn, sOn, tol);
%     sSnsSGn = saturationFraction(sSn, sGn, tol);
%     sSnsNn = saturationFraction(sSn, sOn + sGn, tol);
%     sGnsNn = saturationFraction(sGn, sOn + sSn, tol);
%     
%     sOnsOSn = saturationFraction(sOn, sSn, tol);
%     sGnsSGn = saturationFraction(sGn, sSn, tol);
%     sOnsNn = saturationFraction(sOn, sGn + sSn, tol);


%     muMOS = muO.*muS./(sOnsOSn.*muSa + sSnsOSn.*muOa).^4;
%     muMSG = muS.*muG./(sSnsSGn.*muGa + sGnsSGn.*muSa).^4;
%     muM   = muO.*muS.*muG./(sOnsNn.*muSa.*muGa ...
%                               + sSnsNn.*muOa.*muGa ...
%                               + sGnsNn.*muOa.*muSa).^4;
% 
%     muMOS(isnan(double(muMOS))) = muS(isnan(double(muMOS)));
%     muMSG(isnan(double(muMSG))) = muS(isnan(double(muMSG)));
%     muM(isnan(double(muM))) = muS(isnan(double(muM)));
   
%     muMOS = muO.*muS./((1-sSnsOSn).*muSa + sSnsOSn.*muOa).^4;
%     muMSG = muS.*muG./(sSnsSGn.*muGa + (1-sSnsSGn).*muSa).^4;
%     muM   = muO.*muS.*muG./((1-(sGnsNn + sSnsNn)).*muSa.*muGa ...
%                           + sSnsNn.*muOa.*muGa ...
%                           + sGnsNn.*muOa.*muSa).^4;


%     tol = 1e-5;
%     sGf = saturationFraction(sGn, sOn, tol);
    
%     sSsN_Seff = max((muSmuG.*sGf + muSmuO.*(1-sGf) - (muS./muS_eff).^a)...
%                ./(muSmuG.*sGf + muSmuO.*(1-sGf) - 1),0);

%     sO = max(sO,0);
%     sG = max(sG,0);
%     sS = max(sS,0);

%     tol = 10*eps;
%     sO(sO < tol) = 0;
%     sG(sG < tol) = 0;
%     sS(sS < tol) = 0;
           
                
% %         tol = 1e-2;
% %         tol = eps;
%         sSnsOSn = saturationFraction(sSn, sOn, tol);
%         sSnsSGn = saturationFraction(sSn, sGn, tol);
%         sOnsNn = saturationFraction(sSn, sOn + sGn, tol);
%         sGnsNn = saturationFraction(sGn, sOn + sSn, tol);
% 
%         muMOS = muO.*muS./((1-sSnsOSn).*muSa + sSnsOSn.*muOa).^4;
%         muMSG = muS.*muG./(sSnsSGn.*muGa + (1-sSnsSGn).*muSa).^4;
%         muM   = muO.*muS.*muG./((1-(sGnsNn + sOnsNn)).*muSa.*muGa ...
%                               + sOnsNn.*muOa.*muGa ...
%                               + sGnsNn.*muOa.*muSa).^4;
%     tol = 1e-5;
%     sOsN = saturationFraction(sO, sG + sS, tol);
%     sGsN = saturationFraction(sG, sO + sS, tol);

%         muMOS = muO.*muS.*((sOn + sSn)./(sOn.*muSa + sSn.*muOa)).^4;
%         muMSG = muS.*muG.*((sSn + sGn)./(sSn.*muGa + sGn.*muSa)).^4;
%         muM   = muO.*muG.*muS.*((sOn + sGn + sSn)./...
%                      (sOn.*muSa.*muGa + sGn.*muOa.*muSa + sSn.*muOa.*muGa)).^4;
% 
%         % In the case of sOn = 0 and sSn = 0 etc., we set mixed viscosities to
%         % the unmixed.
%         muMOS(isnan(double(muMOS))) = muO(isnan(double(muMOS)));
%         muMSG(isnan(double(muMSG))) = muG(isnan(double(muMSG)));
%         muM(isnan(double(muM)))     = muS(isnan(double(muM)));
%         sSsN_Seff = max((muSmuG.*sGn + muSmuO.*sOn - (sOn + sGn).*(muS./muS_eff).^a)...
% %                   ./(muSmuG.*sGn + muSmuO.*sOn - (sOn + sGn)),0);
%     if useFrac
%         
%         sGf = saturationFraction(sGn, sOn, tol);
% 
%         sSsN_Seff = max((muSmuG.*sGf + muSmuO.*(1-sGf) - (muS./muS_eff).^a)...
%                    ./(muSmuG.*sGf + muSmuO.*(1-sGf) - 1),0);