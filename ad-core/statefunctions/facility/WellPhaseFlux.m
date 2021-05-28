classdef WellPhaseFlux < StateFunction
    % Get phase-flux between well-bore and reservoir
    properties
    allowCrossFlow = true;
    end

    methods
        function gp = WellPhaseFlux(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'FacilityWellMapping', ...
                               'PressureGradient', ...
                               'WellIndex', ...
                               'Mobility'});
            gp.label = 'q_\alpha';
        end

        function q_ph = evaluateOnDomain(prop, model, state)
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            [dp, wi] = prop.getEvaluatedDependencies(state, 'PressureGradient', 'WellIndex');
            mobw = prop.getEvaluatedDependencies(state, 'Mobility');

            Tdp = -wi.*dp;

            q_ph = calculatePhaseRate(Tdp, mobw, map, prop.allowCrossFlow, model);
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
