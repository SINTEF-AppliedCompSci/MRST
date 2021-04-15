classdef PhasePotentialDifference < StateFunction
    % The potential difference over each interface, for each phase. This is
    % essentially the phase gradient plus any contribution from gravity
    % acting on the averaged phase density.
    properties
        hasGravity; % Is there gravity?
    end
    
    methods
        function gp = PhasePotentialDifference(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('PressureGradient');
            gp.hasGravity = norm(model.getGravityGradient(), inf) > 0;
            if gp.hasGravity
                gp = gp.dependsOn('GravityPotentialDifference');
            end
            gp.label = '\Theta_\alpha';
        end
        function v = evaluateOnDomain(prop, model, state)
            dp = prop.getEvaluatedDependencies(state, 'PressureGradient');
            v = dp;
            if prop.hasGravity
                rhogdz = prop.getEvaluatedDependencies(state, 'GravityPotentialDifference');
                for i = 1:numel(dp)
                    v{i} = v{i} + rhogdz{i};
                end
            end
        end
    end
end

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
