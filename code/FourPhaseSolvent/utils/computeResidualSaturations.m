function [sWres, sOres, sSGres] = computeResidualSaturations(fluid, p, sG, sS)
% Calculate effective residual saturations

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

    % Residual saturations for the immiscible and miscible extrema
    sOres_m    = fluid.sOres_m ;
    sOres_i    = fluid.sOres_i ;
    sSGres_m   = fluid.sSGres_m;
    sSGres_i   = fluid.sSGres_i;
    
    % Misscibility is a function of the solvent fraction in the total gas
    % phase
    M = fluid.Msat(sG, sS).*fluid.Mpres(p);
    
    % Interpolated water/oil residual saturations
    sWres  = fluid.sWres;
    sOres  = M.*sOres_m  + (1 - M).*sOres_i ;
    sSGres = M.*sSGres_m + (1 - M).*sSGres_i;
    
end