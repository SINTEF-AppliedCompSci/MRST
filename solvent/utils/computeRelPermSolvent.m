function [krW_eff, krO_eff, krG_eff, krS_eff] = computeRelPermSolvent(model, p, sW, sO, sG, sS, sWcon, sOr, sGc, mobMult)
% Calulates effective relative permeabilities.

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

    %% Get model fluid

    fluid = model.fluid;
    
    %% Immiscible relative permeabilities
    
    sOr_i   = fluid.sOr_i;
    sGc_i   = fluid.sGc_i;
    srtot   = sWcon + sOr + sGc;
    srtot_i = sWcon + sOr_i + sGc_i;
    
    if model.dynamicEndPointScaling
        sWn  = scaledSaturation(fluid, sW          , srtot, srtot_i, sWcon    , sWcon);
        sOn  = scaledSaturation(fluid, sO          , srtot, srtot_i, sOr      , sOr_i);
        sGTn = scaledSaturation(fluid, sG + sS     , srtot, srtot_i, sGc      , sGc_i);
        sNn  = scaledSaturation(fluid, sO + sG + sS, srtot, srtot_i, sOr + sGc, sOr_i);
    else
        sWn  = sW;
        sOn  = sO;
        sGTn = sG + sS;
        sNn  = sO + sG + sS;
    end

    krW    = fluid.krW(sWn);
    krO_i  = fluid.krOW(sOn); % Oil relperm
    krGT_i = fluid.krG(sGTn); % Total gas relperm

    %% Miscible relative permeabilities
    
    % Mobile saturaitons
    sGTm   = max(sG + sS - sGc, 0);              % Total mobile gas saturation
    sNm    = max(sO + sG + sS - (sOr + sGc), 0); % Total mobile HC saturaion    
    FGTmNm = fluid.satFrac(sGTm, sNm);
    
    % Miscible relperms
    krN    = fluid.krO(sNn);
    krO_m  = fluid.MkrO(FGTmNm).*krN;
    krGT_m = fluid.MkrG(FGTmNm).*krN;

    
    %% Interpolate between the two
    
    % Interpolate between miscible and immiscible cases (water relperm
    % not affected by the solvent)
    M    = fluid.Ms(fluid.satFrac(sS, sG + sS)).*fluid.Mp(p);
    FSGT = fluid.satFrac(sS, sG + sS); % Solvent to total gas saturation fraction
    
    krW_eff  = krW;
    krO_eff  = M.*krO_m  + (1-M).*krO_i;
    krGT_eff = M.*krGT_m + (1-M).*krGT_i;
    krG_eff  = fluid.krFG(1-FSGT).*krGT_eff;
    krS_eff  = fluid.krFS(FSGT  ).*krGT_eff;

    % Modifiy relperm by mobility multiplier (if any)
    krW_eff = mobMult.*krW_eff;
    krO_eff = mobMult.*krO_eff;
    krG_eff = mobMult.*krG_eff;
    krS_eff = mobMult.*krS_eff;
    
    
end

function sn = scaledSaturation(fluid, s, srtot, srtot_i, sr, sr_i)

%     tol = 1e-4;
    sn = max(s-sr,0);
%     sn(sn <= 0) = -tol;
    sn = sn./(1-srtot);
%     sn = min(max(sn,-),1);
    sn = sn.*(1 - srtot_i) + sr_i;
    
end