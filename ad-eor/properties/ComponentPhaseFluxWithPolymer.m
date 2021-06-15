classdef ComponentPhaseFluxWithPolymer < StateFunction
    % Flux of each component, in each phase
    properties (Access = protected)
        mobility_name; % Name of state function where component mobility comes from
    end

    methods
        function cf = ComponentPhaseFluxWithPolymer(model, mob_name)
            cf@StateFunction(model);
            if nargin < 2
                mob_name = 'FaceComponentMobility'; 
            end
            cf.mobility_name = mob_name;
            cf = cf.dependsOn({'PermeabilityPotentialGradient', cf.mobility_name});
            cf.label = 'V_{i,\alpha}';
        end

        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            [kgrad, compMob] = prop.getEvaluatedDependencies(state,...
                'PermeabilityPotentialGradient', prop.mobility_name);
            v = cell(ncomp, nph);
            for c = 1:ncomp
                for ph = 1:nph
                    mob = compMob{c, ph};
                    if ~isempty(mob)
                        v{c, ph} = -mob.*kgrad{ph};
                    end
                end
            end
            

            % when polymer exist, shear effects might play a role
            v = applyShearEffects(v, model, state);
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
