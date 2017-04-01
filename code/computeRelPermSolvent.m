function [krW, krO, krG, krS] = computeRelPermSolvent(sW, sO, sG, sS, p, fluid)

    %% 
    
    krW = @(sW) fluid.krW(sW);
    
    krO_i = fluid.krO(sO);
    krO_m = @(sO, aG, sS) sO./(sO + sG + sS).*fluid.krOW(sO + sG + sS);
    
    
    krGT_i = @(sG, sS) fluid.krG(sG + sS);
    krGT_m = @(sO, sG, sS) (sS + sG)./(sO + sG + sS).*fluid.krOW(sO + sG + sS);
    
    
    krG_i = @(sG, sS) krGT_i(sG, sS).*(sG./(sS + sG));
    krG_m = @(sO, sG, sS) krGT_m(sO, sG, sS).*(sG./(sS + sG));
    
    krS_i = @(sG, sS) krGT_i(sG, sS).*(sS./(sS + sG));
    krS_m = @(sO, sG, sS) krGT_m(sO, sG, sS).*(sS./(sS + sG));

    %%
    
    M = fluid.Msat(sG, sS).*fluid.Mpres(p);
    
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
    
    krW = krW(sW_eff);
    krO = M.*krO_m(sO_eff, sG_eff, sS_eff) + (1-M).*krO_i(sO_eff);
    krG = M.*krG_m(sO_eff, sG_eff, sS_eff) + (1-M).*krG_i(sG_eff, sS_eff);
    krS = M.*krS_m(sO_eff, sG_eff, sS_eff) + (1-M).*krS_i(sG_eff, sS_eff);

end
