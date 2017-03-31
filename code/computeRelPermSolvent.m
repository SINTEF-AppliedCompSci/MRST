function [krW, krO, krG, krS] = computeRelPermSolvent(sW, sO, sG, sS, sS_thres, Nc, fluid)

    %% 
    
    krW = @(sW) fluid.krW(sW);
    
    krO = fluid.krO(sO).*(sS < sS_thres) ...
          + sO./(sO + sG + sS).*fluid.krOW(sO + sG + sS).*(sS >= sS_thres);
    
    krGT = @(sO, sG, sS) fluid.krG(sG + sS).*(sS < sS_thres) ...
         + (sS + sG)./(sO + sG + sS).*fluid.krOW(sO + sG + sS).*(sS >= sS_thres);
    
    krG = @(sO, sG, sS) krGT(sO, sG, sS).*(sG./(sS + sG));
    
    krS = @(sO, sG, sS) krGT(sO, sG, sS).*(sS./(sS + sG));

    %%
    
    
    is_sol = sS > 0;
    M = fluid.Msat(sG, sS);
    
%     M = 0*c;
%     if nnz(is_sol) > 0
%        logNc = log(Nc(is_sol))/log(10);
%        % We cap logNc (as done in Eclipse)
%        logNc = min(max(-20, logNc), 20);
%        M(is_sol) = fluid.miscfact(logNc, 'cellInx', find(is_sol));
%     end
% 
    sOres_m    = fluid.sOres_m ;    % Residual oil saturation     without solvent
    sOres_i    = fluid.sOres_i ;    % Residual oil saturation     without solvent
    sSGres_m   = fluid.sSGres_m;    % Residual oil saturation     without solvent
    sSGres_i   = fluid.sSGres_i;    % Residual oil saturation     without solvent
    sWres = fluid.sWres;
    
    % Interpolated water/oil residual saturations
    sOres  = M.*sOres_m + (1 - M).*sOres_i;
    sSGres = M.*sSGres_m + (1 - M).*sSGres_i;
    
    sW_eff = (sW - sWres )./(1 - sOres - sSGres - sWres );
    sO_eff = (sO - sOres )./(1 - sWres - sSGres - sOres );
    sG_eff = (sG - sSGres)./(1 - sWres - sOres  - sSGres);
    sS_eff = (sS - sSGres)./(1 - sWres - sOres  - sSGres);
    
    sN_eff = sO_eff + sG_eff + sS_eff;
    
    krOW_m = fluid.krOW_m;
    krOW_i = fluid.krOW_i;
    
    krW = krW(sW);
    krO = sO_eff./sN_eff.*(M.*krOW_m(sN_eff) + (1-M).*krOW_i(sN_eff));
    krGT = (sS_eff + sG_eff)./sN_eff.*(M.*krOW_m(sN_eff) + (1-M).*krOW_i(sN_eff));
    krS = krGT.*(sS_eff./(sS_eff + sG_eff));
    krG = krGT.*(sG_eff./(sS_eff + sG_eff));

end
