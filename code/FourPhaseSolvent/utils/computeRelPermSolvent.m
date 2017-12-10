function [krW_eff, krO_eff, krG_eff, krS_eff] = computeRelPermSolvent(fluid, p, sW, sO, sG, sS, sWres, sOres, sSGres, mobMult)
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

    krN    = fluid.krN  (sO + sG + sS, sWres, sOres, sSGres);
    krO_i  = fluid.krO_i(sO          , sWres, sOres, sSGres);
    krGT_i = fluid.krGT_i(sG + sS    , sWres, sOres, sSGres);
              
    M = fluid.Msat(sG, sS).*fluid.Mpres(p);

    sSsGT = sS./(sS + sG);
    sSsGT(isnan(double(sSsGT))) = 0;
    
    % Immiscible gas and solvent relperms
    krG_i = (1-sSsGT).*krGT_i;
    krS_i = sSsGT.*krGT_i;

    
    sOn  = max(sO - sOres,0);
    sGn  = max(sG - sSGres,0);
    sSn  = max(sS - sSGres,0);
    sGTn = max(sG + sS - sSGres, 0);
    sNn  = max(sO + sG + sS - (sOres + sSGres), 0);
    
    sOnsNn = sOn./sNn;
    sOnsNn(isnan(double(sOnsNn))) = 0;
    sGTnsNn = 1 - sOnsNn;
    
    sGnsGTn = sGn./sGTn;
    sGnsGTn(isnan(double(sGnsGTn))) = 0;
    sSnsGTn = 1 - sGnsGTn;
    
    % Miscible relperms
    krO_m = sOnsNn.*krN;
    krGT_m = sGTnsNn.*krN;
    
    krG_m = sGnsGTn.*krGT_m;
    krS_m = sSnsGTn.*krGT_m;
    
    % Interpolate between miscible and immiscible cases (water relperm
    % not affected by the solvent)
    krW_eff = fluid.krW_i(sW, sWres, sOres, sSGres);
    
    krO_eff = M.*krO_m + (1-M).*krO_i;
    krG_eff = M.*krG_m + (1-M).*krG_i;
    krS_eff = M.*krS_m + (1-M).*krS_i;
    
    % Modifiy relperm by mobility multiplier (if any)
    
    % Multiply by sn>0 to avoid erroneous positive relperms due to
    % roundoff-errors in evaluation of miscibility M.
    krW_eff = mobMult.*krW_eff;
    krO_eff = mobMult.*krO_eff.*(sOn>0);
    krG_eff = mobMult.*krG_eff.*(sGn>0);
    krS_eff = mobMult.*krS_eff.*(sSn>0);
    
end