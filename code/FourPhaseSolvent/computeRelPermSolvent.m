function [krW_eff, krO_eff, krG_eff, krS_eff] = computeRelPermSolvent(fluid, p, sW, sO, sG, sS, sWres, sOres, sSGres, mobMult)
    % Calulates effective relative permeabilities.
    

    
%     % Relperm of hydrocarbon to water as a funciton of total HC saturation
%     krN = fluid.krO(sO + sG + sS);
%     
%     % Immiscible relperms for oil and total gas
%     krO_i = fluid.krO(sO);
%     krGT_i = fluid.krG(sG + sS);

    n = 2;

    sr_tot = sWres + sOres + sSGres;
    
    krN = coreyRelperm(sO + sG + sS, ...
                  n, ...
                  sOres + sSGres, ...
                  fluid.krO(1 - sWres), ...
                  sr_tot);
              
    krO_i = coreyRelperm(sO, ...
                  n, ...
                  sOres, ...
                  fluid.krO(1 - (sWres + sSGres)), ...
                  sr_tot);
              
    krGT_i = coreyRelperm(sG + sS, ...
                  n, ...
                  sSGres, ...
                  fluid.krG(1 - (sWres + sOres)), ...
                  sr_tot);


    M = fluid.Msat(sG, sS).*fluid.Mpres(p);
    
%     % Add small value to avoid 0/0-type expressions
%     tol = eps;
%     sO = sO + (sO < tol).*tol;
%     sG = sG + (sG < tol).*tol;
%     sS = sS + (sS < tol).*tol;
%     
%     tol = eps;
%     sGsS = (sG + (sG < tol).*tol)./(sG + (sG < tol).*tol + sS);

    sGsGT = saturationFraction(sG, sG + sS);
    
    % Immiscible gas and solvent relperms
    krG_i = sGsGT.*krGT_i;
    krS_i = (1-sGsGT).*krGT_i;

    
%     % Immiscible gas and solvent relperms
%     krG_i = sG./(sS + sG).*krGT_i;
%     krS_i = sS./(sS + sG).*krGT_i;
    
    sOn = max(sO - sOres, 0);
    sGn = max(sG - sSGres,0);
    sSn = sS;
    sNn = max(sO + sG + sS - (sOres + sSGres), 0);
    sGTn = max(sG + sS - sSGres,0);
    
    sOnsNn  = saturationFraction(sOn, sNn);
    sGTnsNn = saturationFraction(sGTn, sNn);
    sGnsGTn = saturationFraction(sGn, sGTn);
    sSnsGTn = saturationFraction(sSn, sGTn);
    
    % Miscible relperms
    krO_m = sOnsNn.*krN;
    krGT_m = sGTnsNn.*krN;
    krG_m = sGnsGTn.*krGT_m;
    krS_m = sSnsGTn.*krGT_m;
    
%     % Miscible relperms
%     krO_m = sO./(sO + sG + sS).*krN;
%     krGT_m = (sG + sS)./(sO + sG + sS).*krN;
%     krG_m = sG./(sS + sG).*krGT_m;
%     krS_m = sS./(sS + sG).*krGT_m;
    
    % Interpolate between miscible and immiscible cases (water relperm
    % not affected by the solvent)
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

function kr = coreyRelperm(s, n, sr, kwm, sr_tot)
    den = 1 - sr_tot;
    sat = ((s - sr)./den);
    if isa(sat, 'ADI')
        sat.val = max(min(sat.val, 1), 0);
    else
        sat = max(min(sat, 1), 0);
    end
    
    kr = kwm.*sat.^n;
end

function f = saturationFraction(sNum, sDen)

    tol = 10*eps;
    f = sNum./(sDen + (sDen < tol).*tol);
    
end
    