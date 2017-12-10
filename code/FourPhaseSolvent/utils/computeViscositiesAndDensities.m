function [muW_eff , muO_eff , muG_eff , muS_eff , ...
          rhoW_eff, rhoO_eff, rhoG_eff, rhoS_eff, ...
          bW_eff  , bO_eff  , bG_eff  , bS_eff  , ...
          pW      , pG                          ] ...
          = computeViscositiesAndDensities(model, p , sO , sG , sS , sOres , sSGres, rs, rv, isSatO, isSatG)
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

%% Unmixed viscosites and densities at reservoir conditions

    fluid  = model.fluid;
    disgas = isprop(model, 'disgas') && model.disgas;
    vapoil = isprop(model, 'vapoil') && model.disgas;
    
    % Water
    pcOW = 0;
    if isfield(fluid, 'pcOW') && ~isempty(sW)
        pcOW  = fluid.pcOW(sW);
    end
    pW  = p - pcOW;
    bW  = fluid.bW(pW);
    muW = fluid.muW(pW);
    
    rhoW = bW.*fluid.rhoWS;
    
    % Oil    
    if disgas
        bO_i = fluid.bO(p,  rs, isSatO);
        muO  = fluid.muO(p, rs, isSatO);
    else
        bO_i  = fluid.bO(p);
        if isfield(fluid, 'BOxmuO')
            muO = fluid.BOxmuO(p).*bO_i;
        else
            muO = fluid.muO(p);
        end
    end
    
    if any(bO_i < 0)
        warning('Negative oil compressibility present!')
    end
    
    rhoO = bO_i.*(rs*fluid.rhoGS + fluid.rhoOS);
    
    % Gas
    pcOG = 0;
    if isfield(fluid, 'pcOG') && ~isempty(sG)
        Mp   = fluid.Mpres(p);
        pcOG = Mp.*fluid.pcOG(sG) + (1-Mp).*fluid.pcOG(sG + sS);
    end
    pG = p + pcOG;

    if vapoil
        bG_i  = fluid.bG(pG, rv, isSatG);
        muG = fluid.muG(pG, rv, isSatG);
    else
        bG_i = fluid.bG(pG);
        muG = fluid.muG(pG);
    end
    
    if any(bG_i < 0)
        warning('Negative gas compressibility present!')
    end
    
    rhoG = bG_i.*(rv*fluid.rhoOS + fluid.rhoGS);
    
    % Solvent
    bS_i  = fluid.bS(pG);
    muS  = fluid.muS(pG);
    
    rhoS = bS_i.*fluid.rhoSS;
    
    %% Effective viscosites
    
    tol = 1e-8;
%     is_solvent = sS > tol;
    is_solvent = true;

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
    sOnsNn  = sOn./(sOn + sGn + sSn);
    sOnsNn(isnan(double(sOnsNn))) = 0;
    sGnsNn  = sGn./(sOn + sGn + sSn);
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
    muO_eff = muO_eff.*is_solvent + muO.*(~is_solvent);
    muG_eff = muG.^(1-omega).*muMSG.^omega;
    muG_eff = muG_eff.*is_solvent + muG.*(~is_solvent);
    muS_eff = muS.^(1-omega).*muM.^omega;
    muS_eff = muS_eff.*is_solvent + muS.*(~is_solvent);
    
    
    %% Effective densities
    
    % Effective fractional saturations    
    muOmuS = (muO./muS).^a;
    muSmuO = 1./muOmuS;
    muSmuG = (muS./muG).^a;
    sOsN_Oeff = (muOmuS - (muO./muO_eff).^a)./(muOmuS-1);
    sOsN_Oeff(~isfinite(double(sOsN_Oeff))) = 0; 
    sOsN_Geff = (muSmuG - (muS./muG_eff).^a)./(muSmuG-1);
    sOsN_Geff(~isfinite(double(sOsN_Geff))) = 0; 
    
    sGf = sGn./(sOn + sGn).*(sGn>0);
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
    rhoO_eff = rhoO_eff.*is_solvent + rhoO.*(~is_solvent);
                                   
    rhoG_eff = (sOsN_Geff.*rhoS + (1-sOsN_Geff).*rhoG).*(~eq) ...
                                       + ((1-omega)*rhoG + omega*rhoM).*eq;
    rhoG_eff = rhoG_eff.*is_solvent + rhoG.*(~is_solvent);
    
    rhoS_eff = (sSsN_Seff.*rhoS + (1-sSsN_Seff).*(rhoG.*sGf + rhoO.*(1-sGf))).*(~eq) ...
                                       + ((1-omega)*rhoS + omega*rhoM).*eq;
    rhoS_eff = rhoS_eff.*is_solvent + rhoS.*(~is_solvent);
                                   
    % Effective formation volume factors are interpolated using a
    % pressure-dependent miscibility funciton
    
    % Miscible formation volume factors b = rho_eff/rhoS
    bO_m = rhoO_eff./fluid.rhoOS;
    bG_m = rhoG_eff./fluid.rhoGS;
    bS_m = rhoS_eff./fluid.rhoSS;
    
    % Effective formation volume factors are interpolated using a
    % pressure-dependent miscibility funciton
    Mp     = fluid.Mpres(p);
    bW_eff = bW;
    bO_eff = bO_m.*Mp + bO_i.*(1-Mp);
    bG_eff = bG_m.*Mp + bG_i.*(1-Mp);
    bS_eff = bS_m.*Mp + bS_i.*(1-Mp);
    
end