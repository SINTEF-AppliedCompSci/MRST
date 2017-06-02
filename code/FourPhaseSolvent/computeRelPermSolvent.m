function [krW_eff, krO_eff, krG_eff, krS_eff] = computeRelPermSolvent(fluid, p, sW, sO, sG, sS, sWres, sOres, sSGres, mobMult)
% Calulates effective relative permeabilities.

    n = 2;
    sres_tot = sWres + sOres + sSGres;
    
    % Relperm of hydrocarbon to water as a funciton of total HC saturation
    krN = coreyRelperm(sO + sG + sS, ...
                  n, ...
                  sOres + sSGres, ...
                  fluid.krO(1 - sWres), ...
                  sres_tot);

              
    % Immiscible relperms for oil and total gas          
    krO_i = coreyRelperm(sO, ...
                  n, ...
                  sOres, ...
                  fluid.krO(1 - (sWres + sSGres)), ...
                  sres_tot);
    krGT_i = coreyRelperm(sG + sS, ...
                  n, ...
                  sSGres, ...
                  fluid.krG(1 - (sWres + sOres)), ...
                  sres_tot);

    M = fluid.Msat(sG, sS).*fluid.Mpres(p);


    
    sSsGT = saturationFraction(sS, sG);
    
    % Immiscible gas and solvent relperms
    krG_i = (1-sSsGT).*krGT_i;
    krS_i = sSsGT.*krGT_i;

    sOn = max(sO - sOres, 0);
    sGn = max(sG - sSGres,0);
    sSn = max(sS,0);
    sGTn = sGn + sSn;
    
    sOnsNn  = saturationFraction(sOn, sGTn);
    sGTnsNn = 1 - sOnsNn;
    sGnsGTn = saturationFraction(sGn, sSn);
    sSnsGTn = 1 - sGnsGTn;

    % Miscible relperms
    krO_m = sOnsNn.*krN;
    krGT_m = sGTnsNn.*krN;
    krG_m = sGnsGTn.*krGT_m;
    krS_m = sSnsGTn.*krGT_m;
    
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
    
%     sat = max(min((s - sr)./den,1),0);
%     tol = 100*eps;
%     tol = 1e-5;
    tol = 0;
    sat = max(min((s-sr)./(den + tol),1),0);
    
%     if isa(sat, 'ADI')
%         sat.val = max(0,min(sat.val, 1));
%     else
%         sat = max(min(sat, 1), 0);
%     end
    
    kr = kwm.*sat.^n;
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
%     tol = 10*eps;
%     zz = sNum < tol & sDen < tol;
%     sNum = sNum - (sNum < tol).*sNum;
%     sDen = sDen - (sDen < tol).*sDen;
%     f = sNum./(sDen + (sDen < tol).*tol).*(~zz) + 0.*zz;
%        
% end
    