function [vW, mobW] = getFluxAndPropsWater(model, pW, krW, T, gdz)
%  Function to compute the water fluxes and mobilities.
% 
% This function is modified from a file in The MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/ad-eor/utils/getFluxAndPropsWaterPolymer_BO.m 
%
% We refer to that function for a complete commented version of the file. 

%{
Partial copyright 2009-2020, SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2020, NORCE Norwegian Research Centre AS, Computational 
Geosciences and Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}
    fluid = model.fluid;
    s = model.operators;
    dpW    = s.Grad(pW) - fluid.rhoWS*gdz;   
    upcw  = (value(dpW)<=0);
    [krWf, krW] = s.splitFaceCellValue(s, upcw, krW);    
    mobW   = krW/fluid.muW;    
    vW = -(krWf/fluid.muW).*T.*dpW;
end