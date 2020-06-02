function state = addDerivedQuantitiesDC(model, state)
%
% SYNOPSIS:
%   function state = addDerivedQuantitiesDC(model, state)
%
% DESCRIPTION: 
%   Computes extra mechanical volumetric strain fields for the DP model 
%   (intrinsic matrix and fracture strains) from the primary variable of 
%   the mechanical state (state.xd) and the pressure variables (state.pressure
%   and state.pressure_matrix). 
%
% PARAMETERS:
%   model - Poromechanical model
%   state - State variable
%
% RETURNS:
%   state - State variable
%
% EXAMPLE:
%
% SEE ALSO:
%
%{
Copyright 2009-2020 SINTEF ICT, Applied Mathematics.

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
    % Extract quantities
    p_f = model.fluidModel.getProp(state, 'pressure');
    p_m = model.fluidModel.getProp(state, 'pressure_matrix');
    total_vol = prod(max(model.G.nodes.coords));
    strain = model.mechModel.getProp(state, 'strain'); % cell-wise strain
    global_strain = strain.*(model.G.cells.volumes./total_vol);
    
    % Define parameters required for calculation
    cM = model.mechModel.constitutive_coefficients_object;  
    v_m = model.rock_matrix.vol_fraction; % fracture volume fraction
    invb_m = 1./cM.iP.b_m(:,1);                   
    
    % Calculate intrinsic volumetric strains. Calculation based on comparing expressions
    % for effective change in Lagrangian matrix porosity and the volumetric
    % average of the local change of Lagrangian matrix porosity.
    v_part = model.G.cells.volumes./total_vol; % partition volume
    vol_strain_mat = (invb_m).*((1/v_m)*(sum(cM.B_m.*global_strain, 2) + cM.invN_m.*p_m.*v_part ...
                                + cM.invQ.*p_f.*v_part) - cM.iP.invn_m.*p_m.*v_part);
    vol_strain_frac = (model.mechModel.operators.trace(global_strain)-v_m.*vol_strain_mat)./(1-v_m);
    
    state = setProp(model.mechModel, state, 'vol_strain_frac', vol_strain_frac);
    state = setProp(model.mechModel, state, 'vol_strain_mat', vol_strain_mat);      
end
