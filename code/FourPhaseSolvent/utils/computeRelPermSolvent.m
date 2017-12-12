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
    
    krO_i  = fluid.krO(sO);      % Oil relperm
    krGT_i = fluid.krG(sG + sS); % Total gas relperm
    
    sS(sS < fluid.smin) = 0;
    FSG = fluid.satFrac(sS, sG + sS); % Solvent to total gas saturation fraction
    krG_i = fluid.krFG((1-FSG)).*krGT_i;
    krS_i = fluid.krFS(FSG    ).*krGT_i;

    %% Miscible relative permeabilities
    
    % Mobile saturaitons
    sOn  = max(sO - sOr,0);
    sGn  = max(sG - sGc,0);
    sGTn = max(sG + sS - sGc, 0);              % Total mobile gas saturation
    sNn  = max(sO + sG + sS - (sOr + sGc), 0); % Total mobile HC saturaion
    
    FOnNn = fluid.satFrac(sOn, sNn);   % Oil to HC saturation fraction
    FGnGTn = fluid.satFrac(sGn, sGTn); % Gas to HC saturation fraction
    
    % Miscible relperms
    krN    = fluid.krOW(sO + sG + sS);
    krO_m  = fluid.MkrO(FOnNn).*krN;
    krGT_m = fluid.MkrG(1-FOnNn).*krN;
    krG_m  = fluid.krFG(FGnGTn).*krGT_m;
    krS_m  = fluid.krFS(1-FGnGTn).*krGT_m;
    
    %% Interpolate between the two
    
    % Interpolate between miscible and immiscible cases (water relperm
    % not affected by the solvent)
    M = fluid.Ms(fluid.satFrac(sS, sG + sS)).*fluid.Mp(p);
    
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