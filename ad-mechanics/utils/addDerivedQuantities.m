function state = addDerivedQuantities(model, state)
%
%
% SYNOPSIS:
%   function state = addDerivedQuantities(model, state)
%
% DESCRIPTION:
%   Computes extra mechanical fields (such as strain, stress) from the
%   primary variable of the mechanical state (state.xd)
%
% PARAMETERS:
%   model - Mechanical model
%   state - State variable
%
% RETURNS:
%   state - State variable

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

    u = model.operators.V_dir; % zeros(G.gridim * G.nodes.num, 1);
    u(~model.operators.isdirdofs) = getProp(model, state, 'xd');
    state = setProp(model, state, 'u', u);
    uu = reshape(u, model.G.griddim, [])';
    state = setProp(model, state, 'uu', uu);
    
    % calculate div, stress, eigen values + + +
    vdiv = model.operators.ovol_div * u ./ model.G.cells.volumes;
    state = setProp(model, state, 'vdiv', vdiv);
    
    if(model.G.griddim == 2)
        lin_dim = 3; % dimension of the non trivial linear degree of freedom which
    elseif(model.G.griddim == 3)
        lin_dim = 6;
    else
        error('Wrong dimsension');
    end
    stress = reshape(model.operators.global_stress * u, lin_dim, [])';
    strain = reshape(model.operators.global_strain * u, lin_dim, [])';
    
    state = setProp(model, state, 'stress', stress);
    state = setProp(model, state, 'strain', strain);
    
    
end
