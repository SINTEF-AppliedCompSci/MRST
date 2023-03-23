classdef MpfaPhaseUpwindFlag < StateFunction
    % Phase-upstream flag for each phase
    properties

    end
    
    methods
        function gp = MpfaPhaseUpwindFlag(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('PermeabilityPotentialGradient');
            gp.label = '\Theta_\alpha \leq 0';
        end

        function flags = evaluateOnDomain(prop, model, state)
            v = model.getProp(state, 'PermeabilityPotentialGradient');
            nph = numel(v);
            flags = cell(1, nph);
            for i = 1:nph
                flags{i} = v{i} <= 0;
            end
        end
    end
end

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}

