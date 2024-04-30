classdef CO2VEDissolvedFlux < StateFunction 

    properties (Access = protected)
    end

    methods
        function fm = CO2VEDissolvedFlux(model)
            fm@StateFunction(model);
            fm = fm.dependsOn('ComponentPhaseFlux');
            fm.label = 'v_{rs}';
        end

        function v = evaluateOnDomain(prop, model, state)
            c_phase_fluxes = prop.getEvaluatedDependencies(state, 'ComponentPhaseFlux');
            
            water_phase_ix = model.getPhaseIndex('W');
            co2_phase_ix = model.getPhaseIndex('G');
            
            water_component_ix = 1; % @@ hard-coded for now
            co2_component_ix = 2;   % @@
            
            v = c_phase_fluxes{co2_component_ix, water_phase_ix};
        end
    end
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
