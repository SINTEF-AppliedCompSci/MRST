classdef PhaseFlux < StateFunction
    % Phase flux for each phase. The volumetric flux associated with each
    % internal interface.
    properties

    end
    
    methods
        function fm = PhaseFlux(model)
            fm@StateFunction(model);
            fm = fm.dependsOn({'FaceMobility', 'PermeabilityPotentialGradient'});
            fm.label = 'v_\alpha';
        end

        
        function v = evaluateOnDomain(prop, model, state)
            [mob, kgrad] = prop.getEvaluatedDependencies(state,...
                'FaceMobility', 'PermeabilityPotentialGradient');
            nph = numel(mob);
            v = cell(1, nph);
            for i = 1:nph
                v{i} = -mob{i}.*kgrad{i};
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
