function [krW_eff, krO_eff, krG_eff, krS_eff] = computeRelPermSolvent(model, p, sW, sO, sG, sS, sWr, sOr, sGc, mobMult)
% Calulates effective relative permeabilities.

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

    %% Get model fluid

    fluid = model.fluid;
    
    %% Immiscible relative permeabilities
    
%     sS(sS < fluid.smin) = 0;
%     sG(sG < fluid.smin) = 0;
    
    sOr_i = fluid.sOr_i; sGc_i = fluid.sGc_i;
    srtot = sWr + sOr + sGc;
    srtot_i = sWr + sOr_i + sGc_i;
    
    sWn  = scaledSaturation(fluid, sW, srtot, srtot_i, sWr, sWr  );
    sOn  = scaledSaturation(fluid, sO, srtot, srtot_i, sOr, sOr_i);
    sGTn = scaledSaturation(fluid, sG + sS, srtot, srtot_i, sGc, sGc_i);
    
%     sOn = scaledSaturation(fluid, 1-(sW + sG), srtot, srtot_i, sOr, sOr_i);
%     sWn = sW;
%     sOn = sO;
%     sGTn = sG + sS;
    
    krW    = fluid.krW(sWn);
    krO_i  = fluid.krO(sOn);  % Oil relperm
    krGT_i = fluid.krG(sGTn); % Total gas relperm
    
    
    FSG   = fluid.satFrac(sS, sG + sS); % Solvent to total gas saturation fraction
    krG_i = fluid.krFG((1-FSG)).*krGT_i;
    krS_i = fluid.krFS(FSG    ).*krGT_i;

    %% Miscible relative permeabilities
    
    % Mobile saturaitons
    sOm  = max(sO - sOr, 0);
    sGm  = max(sG - sGc, 0);
    sGTm = max(sG + sS - sGc, 0);              % Total mobile gas saturation
    sNm  = max(sO + sG + sS - (sOr + sGc), 0); % Total mobile HC saturaion
    
    FOmNm = fluid.satFrac(sOm, sNm);   % Oil to HC saturation fraction
    FGmGTm = fluid.satFrac(sGm, sGTm); % Gas to HC saturation fraction
    
    FGTmNm = fluid.satFrac(sGTm, sNm);
    
    % Miscible relperms
    sNn = scaledSaturation(fluid, sO + sG + sS, srtot, srtot_i, sOr + sGc, sOr_i);
    
%     sNn = sO + sG + sS;
    krN    = fluid.krOW(sNn);
%     krO_m  = fluid.MkrO(1-FOmNm).*krN;
%     krGT_m = fluid.MkrG(1-FOmNm).*krN;
    krO_m  = fluid.MkrO(FGTmNm).*krN;
    krGT_m = fluid.MkrG(FGTmNm).*krN;
%     krG_m  = fluid.krFG(FGmGTm).*krGT_m;
%     krS_m  = fluid.krFS(1-FGmGTm).*krGT_m;
    krG_m  = fluid.krFG(1-FSG).*krGT_m;
    krS_m  = fluid.krFS(FSG).*krGT_m;
    
    %% Interpolate between the two
    
    % Interpolate between miscible and immiscible cases (water relperm
    % not affected by the solvent)
    M = fluid.Ms(fluid.satFrac(sS, sG + sS)).*fluid.Mp(p);
    
    krW_eff = krW;
    krO_eff = M.*krO_m + (1-M).*krO_i;
    krG_eff = M.*krG_m + (1-M).*krG_i;
    krS_eff = M.*krS_m + (1-M).*krS_i;
    
    % Modifiy relperm by mobility multiplier (if any)
    krW_eff = mobMult.*krW_eff;
    krO_eff = mobMult.*krO_eff;
    krG_eff = mobMult.*krG_eff;
    krS_eff = mobMult.*krS_eff;
    
end

function sn = scaledSaturation(fluid, s, srtot, srtot_i, sr, sr_i)

    sn = max(s - sr, 0);
    sn = fluid.satFrac(sn, 1 - srtot);
    sn = sn.*(1 - srtot_i) + sr_i;
    sn = min(sn, 1);
    
end