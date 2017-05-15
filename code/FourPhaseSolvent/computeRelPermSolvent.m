function [krW_eff, krO_eff, krG_eff, krS_eff] = computeRelPermSolvent(fluid, p, sW, sO, sG, sS, sWres, sOres, sSGres, mobMult)
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

end