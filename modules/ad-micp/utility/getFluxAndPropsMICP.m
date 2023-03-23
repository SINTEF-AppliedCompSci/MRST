function [vW, vM, vO, vU, model] = ...
                         getFluxAndPropsMICP(model, pW, m, o, u, gdz, poro)
% Function to compute the water and component fluxes and water mobility
% 
% This function is modified from a file in The MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/ad-eor/utils/getFluxAndPropsWaterPolymer_BO.m 
%
% We refer to that function for a complete commented version of the file. 

%{
Partial copyright 2009-2021, SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2021, NORCE Norwegian Research Centre AS, Computational 
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
    % Compute current permeability
    K = model.fluid.K(poro);        
    
    % Update rock and operators witn the current permebaility
    model.rock      = makeRock(model.G, K, poro);
    model.operators = setupOperatorsTPFA(model.G, model.rock);
    
    fluid = model.fluid;
    s = model.operators;

    dpW    = s.Grad(pW) - fluid.rhoWS * gdz;    
    upcw  = (value(dpW) <= 0);   
    % Water
    vW = -(1 / fluid.muw) * s.T .* dpW;
    % Microbes
    [mf, ~] = s.splitFaceCellValue(s, upcw, m);
    vM = -(1 / fluid.muw) * mf .* s.T .* dpW;
    % Oxygen
    [of, ~] = s.splitFaceCellValue(s, upcw, o);
    vO = -(1 / fluid.muw) * of .* s.T .* dpW;
    % Urea
    [uf, ~] = s.splitFaceCellValue(s, upcw, u);
    vU = -(1 / fluid.muw) * uf .* s.T .* dpW;
end