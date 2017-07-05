function [muW_eff, muO_eff, muG_eff, muS_eff, rhoW_eff, rhoO_eff, rhoG_eff, rhoS_eff] = computeViscositiesAndDensities(fluid, p , sO , sG , sS , sOres , sSGres )
% Calculates effective viscosities and densities using Todd-Longstaff
% model + 1/4th-power mixing rule

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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
    sOsN_Oeff = (muOmuS - (muO./muO_eff).^a)./(muOmuS-1);
    sOsN_Oeff(~isfinite(double(sOsN_Oeff))) = 0; 
    sOsN_Geff = (muSmuG - (muS./muG_eff).^a)./(muSmuG-1);
    sOsN_Geff(~isfinite(double(sOsN_Geff))) = 0; 
    
    sGf = sGn./(sOn + sGn);
    sGf(isnan(double(sGf))) = 0;

    sSsN_Seff = (muSmuG.*sGf + muSmuO.*(1-sGf) - (muS./muS_eff).^a)...
               ./(muSmuG.*sGf + muSmuO.*(1-sGf) - 1);
    sSsN_Seff(~isfinite(double(sSsN_Seff))) = 0; 
           
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
    
    rhoS_eff = (sSsN_Seff.*rhoS + (1-sSsN_Seff).*(rhoG.*sGf + rhoO.*(1-sGf))).*(~eq) ...
                                       + ((1-omega)*rhoS + omega*rhoM).*eq;
    
end