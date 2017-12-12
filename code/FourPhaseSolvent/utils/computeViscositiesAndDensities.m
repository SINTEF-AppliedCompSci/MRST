function [muW_eff, muO_eff, muG_eff, muS_eff, rhoW_eff, rhoO_eff, rhoG_eff, rhoS_eff, bW_eff, bO_eff, bG_eff, bS_eff, pW, pG] = computeViscositiesAndDensities(model, p , sO , sG , sS , sOr , sGc, rs, rv, isSatO, isSatG)
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
    
    % Caluculate mobile saturations
    sOn = max(sO - sOr, 0);
    sGn = max(sG - sGc, 0);
    sSn = max(sS - sGc, 0);
    
    % Calculate mixed viscosities
    a    = 1/4;
    muSa = muS.^a;
    muGa = muG.^a;
    muOa = muO.^a;

    % Saturation fractions
    FSnOSn = fluid.satFrac(sSn, sOn + sSn);
    FSnGSn = fluid.satFrac(sSn, sGn + sSn);
    FOnNn  = fluid.satFrac(sOn, sOn + sGn  + sSn);
    FGnNn  = fluid.satFrac(sGn, sOn + sGn  + sSn);
    
    % Mixed viscosities
    muMOS = muO.*muS./((1-FSnOSn).*muSa + FSnOSn.*muOa).^(1/a); % Oil-solvent
    muMSG = muS.*muG./(FSnGSn.*muGa + (1-FSnGSn).*muSa).^(1/a); % Solvent-gas
    muM   = muO.*muS.*muG./(FOnNn.*muSa.*muGa           ... 
                          + (1-FOnNn-FGnNn).*muOa.*muGa ...
                          + FGnNn.*muOa.*muSa).^(1/a);          % Oil-gas-solvent
    
    % Effective viscosities are determined by the mixing parameter
    omega = fluid.mixPar;
    muW_eff = muW;
    muO_eff = muO.^(1-omega).*muMOS.^omega;
    muG_eff = muG.^(1-omega).*muMSG.^omega;    
    muS_eff = muS.^(1-omega).*muM.^omega;
    
    %% Effective densities
    
    % Effective fractional saturations    
    muOmuS = (muO./muS).^a;
    muSmuO = 1./muOmuS;
    muSmuG = (muS./muG).^a;
    
    % Expressions are sinuglar if muO == muG, in which case we replace the
    % by a simple interpolation rho*(1-omega) + rhoM*omega
    eq = abs(muO - muS) < 1e-10 | abs(muS - muG) < 1e-10;

    FON_Oeff = (muOmuS - (muO./muO_eff).^a)./(muOmuS-1);
    FON_Oeff(~isfinite(double(FON_Oeff))) = 0; 
    FON_Geff = (muSmuG - (muS./muG_eff).^a)./(muSmuG-1);
    FON_Geff(~isfinite(double(FON_Geff))) = 0; 
    FGnOGn = fluid.satFrac(sGn, sOn + sGn);
    
    FSN_Seff = (muSmuG.*FGnOGn + muSmuO.*(1-FGnOGn) - (muS./muS_eff).^a)...
               ./(muSmuG.*FGnOGn + muSmuO.*(1-FGnOGn) - 1);
    FSN_Seff(~isfinite(double(FSN_Seff))) = 0; 
           
    FON = fluid.satFrac(sO, sO + sG + sS);
    FGN = fluid.satFrac(sG, sO + sG + sS);
    
    rhoM = rhoO.*FON + rhoG.*FGN + rhoS.*(1 - (FON + FGN));
    
    % Calulcate mixed densities
    rhoW_eff = rhoW;
    
    rhoO_eff = (FON_Oeff.*rhoO + (1-FON_Oeff).*rhoS).*(~eq) ...
                                       + ((1-omega)*rhoO + omega*rhoM).*eq;
    rhoG_eff = (FON_Geff.*rhoS + (1-FON_Geff).*rhoG).*(~eq) ...
                                       + ((1-omega)*rhoG + omega*rhoM).*eq;
    rhoS_eff = (FSN_Seff.*rhoS + (1-FSN_Seff).*(rhoG.*FGnOGn + rhoO.*(1-FGnOGn))).*(~eq) ...
                                       + ((1-omega)*rhoS + omega*rhoM).*eq;
                                   
    % Effective formation volume factors are interpolated using a
    % pressure-dependent miscibility funciton
    
    % Miscible formation volume factors b = rho_eff/rhoS
    bO_m = rhoO_eff./fluid.rhoOS;
    bG_m = rhoG_eff./fluid.rhoGS;
    bS_m = rhoS_eff./fluid.rhoSS;
    
    % Effective formation volume factors are interpolated using a
    % pressure-dependent miscibility funciton
    Mp     = fluid.Mp(p);
    bW_eff = bW;
    bO_eff = bO_m.*Mp + bO_i.*(1-Mp);
    bG_eff = bG_m.*Mp + bG_i.*(1-Mp);
    bS_eff = bS_m.*Mp + bS_i.*(1-Mp);
    
end