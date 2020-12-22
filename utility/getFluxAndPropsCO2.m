function [vO, mobO] = getFluxAndPropsCO2(model, pO, krO, T, gdz, varargin)
%  Function to compute the co2 fluxes and mobilities. 
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
    dpO    = s.Grad(pO) - fluid.rhoOS.*gdz;    
    upco  = (value(dpO)<=0);
    [krOf, krO] = s.splitFaceCellValue(s, upco, krO);    
    mobO   = krO/fluid.muO;    
    vO = -(krOf/fluid.muO).*T.*dpO;
end