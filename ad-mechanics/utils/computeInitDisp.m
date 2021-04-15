function state = computeInitDisp(model, state, uu, varargin)
%
%
% SYNOPSIS:
%   function state = computeInitDisp(model, state, uu, varargin)
%
% DESCRIPTION:
%   Compute the initial displacement
%
% PARAMETERS:
%   model    - Model that contains a mechanical model (model.mechModel)
%   state    - Initial state, to be modified
%   uu       - Given initial displacement, if the initial state is not
%              computed by solving the mechanical system
%
% OPTIONAL PARAMETERS:
%  'pressure' - fluid pressure. Used as input for mechanical system.
%
% RETURNS:
%   state - State updated with computed mechanical variables.

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

    opt = struct('mech_equil', true, ...
                 'pressure'  , []   );
    opt = merge_options(opt, varargin{:});

    mechModel = model.mechModel;
    
    if (opt.mech_equil)

        mechmodel = MechanicModel(mechModel.G, mechModel.rock, mechModel.mech);

        state = mechmodel.setProp(state, 'xd', mechmodel.operators.rhs); % Dummy values, just used
                                                                         % to get the correct dimension.
        
        if ~isempty(opt.pressure)
            drivingForces.fluidp = opt.pressure;
        else
            drivingForces.fluidp = 0;
        end

        solver = NonLinearSolver(); % The problem is linear but we use the general
                                    % framework. The cost should be
                                    % minimal.
        [state, failure, report] = ...
            solveMinistep(solver, mechmodel, state, state, 0, drivingForces);

    else
        error('not checked')
        u = reshape(uu', [], 1);
        state.xd = zeros(size(model.operators.A, 2), 1);
        state.xd = u(~model.operators.isdirdofs);
    end

    state = addDerivedQuantities(mechmodel, state);

end
