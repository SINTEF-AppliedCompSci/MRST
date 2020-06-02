function [eqs, names, types] = equationsDCPoroMechanics(x, model, fluidp, fluidp_matrix)
%
% SYNOPSIS:
%   function [eqs, names, types] = equationsDPPoroMechanics(x, model, fluidp, fluidp_matrix)
%
% DESCRIPTION: 
%   Assembles the residual for dual-permeability momentum balance. The
%   function takes fluid input, given by the pore pressures.
%
% PARAMETERS:
%   x                 - Displacement
%   model             - Model class instance that is used
%   fluidp            - Fluid pressure
%   fluidp_matrix     - Fluid pressure of the matrix
%
% RETURNS:
%   eqs   - The residual values as ADI variables (that is with the Jacobian)
%           if the inputs were also ADI.
%   names - The name of each equations
%   types - The type of each equations
%
% EXAMPLE: 
%
% SEE ALSO: DualContinuumMechWaterModel
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

    
    G = model.G;
    d = G.griddim;
    s = model.operators;
    cM = model.constitutive_coefficients_object;
    B_f = repmat(cM.B_f(:,1:d), 1, size(s.ovol_div,2)/d);
    B_f = B_f(:, ~s.isdirdofs);
    B_m = repmat(cM.B_m(:,1:d), 1, size(s.ovol_div,2)/d);
    B_m = B_m(:, ~s.isdirdofs);

    eqs{1} = s.A * x - (s.gradP .* B_m') * fluidp_matrix - (s.gradP .*  B_f') * fluidp - s.rhs;

    % normalization constant
    fac =  1 / (1e6 * mean(G.cells.volumes)); 
    eqs{1} = eqs{1} * fac;
    names = {'disp'};
    types = {'disp_dofs'};

end