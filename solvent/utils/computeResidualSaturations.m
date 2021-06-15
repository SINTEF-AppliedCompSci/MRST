function [sWcon, sOr, sGc] = computeResidualSaturations(model, p, sW, sG, sS, state0)
% Calculate effective residual saturations

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

    fluid = model.fluid;
    
    % Residual saturations for the immiscible and miscible extrema
    sOr_m = fluid.sOr_m(sW);
    sOr_i = fluid.sOr_i;
    sGc_m = fluid.sGc_m(sW);
    sGc_i = fluid.sGc_i;
    
    % Misscibility is a function of the solvent fraction in the total gas
    % phase
    
    FSGT = fluid.satFrac(sS, sG + sS);
    M = fluid.Ms(FSGT).*fluid.Mp(p);
    
    % Interpolated water/oil residual saturations
    sWcon = fluid.sWcon;
    sOr = M.*sOr_m + (1 - M).*sOr_i ;
    sGc = M.*sGc_m + (1 - M).*sGc_i;
        
    if nargin > 5 && isfield(state0, 'sOr') && isfield(state0, 'sGc') && model.hystereticResSat
        sOr = min(sOr, state0.sOr);
        sGc = min(sGc, state0.sGc);
    end

end