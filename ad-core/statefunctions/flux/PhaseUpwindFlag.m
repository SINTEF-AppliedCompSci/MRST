classdef PhaseUpwindFlag < StateFunction
    % Phase-upstream flag for each phase
    properties

    end
    
    methods
        function gp = PhaseUpwindFlag(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('PhasePotentialDifference');
            gp.label = '\Theta_\alpha \leq 0';
        end

        function flags = evaluateOnDomain(prop, model, state)
            pot = prop.getEvaluatedDependencies(state, 'PhasePotentialDifference');
            nph = numel(pot);
            flags = cell(1, nph);
            for i = 1:nph
                flags{i} = pot{i} <= 0;
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
