classdef PermeabilityPotentialGradient < StateFunction
    % The difference in potential over a face, multiplied with a
    % discretized permeability
    properties
        PermeabilityGradientDiscretization
    end
    
    methods
        function pp = PermeabilityPotentialGradient(model, kgrad)
            pp@StateFunction(model);
            pp.PermeabilityGradientDiscretization = kgrad;
            pp = pp.dependsOn('PhasePotentialDifference');
            if isa(kgrad, 'TwoPointFluxApproximation')
                pp = pp.dependsOn('Transmissibility');
            end
            pp.label = 'T_f \Theta_\alpha';
        end
        
        function v = evaluateOnDomain(prop, model, state)
            pot = prop.getEvaluatedDependencies(state, 'PhasePotentialDifference');
            nph = numel(pot);
            kgrad = prop.PermeabilityGradientDiscretization;
            v = cell(1, nph);
            for i = 1:nph
                v{i} = kgrad.getPermeabilityGradient(model, state, pot{i});
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
