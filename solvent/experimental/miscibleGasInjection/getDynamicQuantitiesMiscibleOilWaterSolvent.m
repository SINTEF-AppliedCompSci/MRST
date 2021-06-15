function [kr_eff, mu_eff, rho_eff, b_eff, b0_eff, pvMult, pvMult0, T] ...
  = getDynamicQuantitiesMiscibleOilWaterSolvent(model, p0, p, sW, sO, sG, sO0, sG0)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    fluid = model.fluid;
    op    = model.operators;
    
    % Get multipliers
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);
    
    %% Transmissibility
    T = op.T.*transMult;
    
    %% Relative permeabilites
    
    % Residual saturations
    sOres = fluid.sOres;
    sGres = fluid.sGres;
    
    sOn = max(sO - sOres,0);
    sGn = max(sG - sGres,0);
    
    sOn0 = max(sO0 - sOres,0);
    sGn0 = max(sG0 - sGres,0);
    sNn0 = sOn0 + sGn0;
    
    % Effective relative permeabilites
    krN = fluid.krO(sO + sG);
    krW_eff = fluid.krW(sW);

    % Add smalll value to avoid 0/0-type expressions
    
    tol = eps;
    sOn = sOn + tol*(abs(sOn) < tol);
    sGn = sGn + tol*(abs(sGn) < tol);
    sNn = sOn + sGn;
    
    krO_eff = sOn./sNn.*krN;
    krG_eff = sGn./sNn.*krN;
   
    % Multiply by mobility multiplyer
    krW_eff = mobMult.*krW_eff;
    krO_eff = mobMult.*krO_eff;
    krG_eff = mobMult.*krG_eff;
    
    kr_eff = {krW_eff, krO_eff, krG_eff};
        
    %% Viscosities

    mu_eff  = calculateViscosities(fluid, p , sOn , sGn , sNn );
    mu0_eff = calculateViscosities(fluid, p0, sOn0, sGn0, sNn0);
    
    %% Densities

    rho_eff  = calculateDensities(fluid, p , mu_eff , sOn, sGn, sNn); 
    rho0_eff = calculateDensities(fluid, p0, mu0_eff, sOn, sGn, sNn); 
        
    %% Formation volume factors
    
    b_eff  = calculateFormationVolumeFactors(fluid, p , rho_eff );
    b0_eff = calculateFormationVolumeFactors(fluid, p0, rho0_eff);
        
end 

function mu_eff = calculateViscosities(fluid, p, sOn, sGn, sNn)

 % Unmixed phase viscosities (water not affected by solvent)
    muW = fluid.muW(p);
    muO = fluid.muO(p);
    muG = fluid.muG(p);
    
    % Mixing parameter
    omega = fluid.mixPar;
    
    % Mixed viscosity
    a = 1/4;
    muM = muO.*muG./(sGn./sNn.*muO.^a + sOn./sNn.*muG.^a).^4;
    
    % Effective viscosities
    muW_eff = muW;
    muO_eff = muO.^(1-omega).*muM.^omega;
    muG_eff = muG.^(1-omega).*muM.^omega;
    
    mu_eff = {muW_eff, muO_eff, muG_eff};
    
end

function rho_eff = calculateDensities(fluid, p, mu_eff, sO, sG, sN)
    
    muO = fluid.muO(p);
    muG = fluid.muG(p);
    
    muO_eff = mu_eff{2};
    muG_eff = mu_eff{3};
    
        
    omega = fluid.mixPar;
    
    a = 1/4;

    % Effective saturation fractions
    r = muO./muG;
    sR_Oeff = (r.^a - (muO./muO_eff).^a)./(r.^a - 1);
    sR_Geff = (r.^a - (muO./muG_eff).^a)./(r.^a - 1);

    
    % Unmixed phase viscosities (water not affected by solvent)
    rhoO = fluid.bO(p).*fluid.rhoOS;
    rhoG = fluid.bG(p).*fluid.rhoGS;
    
    rhoM = rhoO.*(sO./sN) + rhoG.*(sG./sN);
    
    % Expressions are sinuglar if muO == muG
    tol = 1e-10;
    eq = abs(muO - muG) < tol;
    
    rhoW_eff = fluid.bW(p).*fluid.rhoWS; 
    rhoO_eff = (rhoO.*sR_Oeff + rhoG.*(1 - sR_Oeff)).*(~eq) ...
                                     + ((1-omega).*rhoO + omega.*rhoM).*eq;
    rhoG_eff = (rhoO.*sR_Geff + rhoG.*(1 - sR_Geff)).*(~eq) ...
                                     + ((1-omega).*rhoG + omega.*rhoM).*eq;
    
    rho_eff = {rhoW_eff, rhoO_eff, rhoG_eff};
    
end

function b_eff = calculateFormationVolumeFactors(fluid, p, rho_eff)

    rhoO_eff = rho_eff{2};
    rhoG_eff = rho_eff{3};

    bW_eff = fluid.bW(p);
    bO_eff = rhoO_eff./fluid.rhoOS;
    bG_eff = rhoG_eff./fluid.rhoGS;
    
    b_eff = {bW_eff, bO_eff, bG_eff};
    
end
