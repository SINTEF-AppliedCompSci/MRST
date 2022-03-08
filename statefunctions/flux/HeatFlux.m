classdef HeatFlux < StateFunction
%State function for computing combined heat flux from conductive and
%advective contributions
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function flux = HeatFlux(model)
            flux@StateFunction(model);
            flux = flux.dependsOn({'ConductiveHeatFlux', 'AdvectiveHeatFlux'});
            flux.label = 'H';
        end
        
        %-----------------------------------------------------------------%
        function v = evaluateOnDomain(prop, model, state)
            [vCond, vAdv] = prop.getEvaluatedDependencies(state,      ...
                                                'ConductiveHeatFlux', ...
                                                'AdvectiveHeatFlux' );
            % Add up fluxes
            nph = model.getNumberOfPhases();
            v   = vCond;
            for i = 1:nph
                v = v + vAdv{i};
            end
        end 
    end
    
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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