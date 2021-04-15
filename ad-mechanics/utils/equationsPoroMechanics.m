function [eqs, names, types] = equationsPoroMechanics(x, model, fluidp)
%
%
% SYNOPSIS:
%   function [eqs, names, types] = equationsPoroMechanics(x, model, fluidp)
%
% DESCRIPTION:
%   Assemble the residual equations for the mechanical system. The
%   function takes fluid input, given by the pore pressure.
%
% PARAMETERS:
%   x      - Displacement
%   model  - Model class instance that is used
%   fluidp - Fluid pressure
%
% RETURNS:
%   eqs   - The residual values as ADI variables (that is with the Jacobian)
%           if the inputs were also ADI.
%   names - The name of each equations
%   types - The type of each equations
%
% SEE ALSO:
%   MechBlackOilModel, MechOilWaterModel, MechWaterModel

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

    G = model.G;
    s = model.operators;
    alpha = model.rock.alpha;

    eqs{1} = s.A * x - s.gradP * (alpha .* fluidp) - s.rhs;

    % normalization constant
    fac =  1 / (1e6 * mean(G.cells.volumes));
    eqs{1} = eqs{1} * fac;
    names = {'disp'};
    types = {'disp_dofs'};

end