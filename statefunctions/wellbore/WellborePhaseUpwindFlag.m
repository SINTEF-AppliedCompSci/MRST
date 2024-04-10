classdef WellborePhaseUpwindFlag < StateFunction
% State function for upwind flag in WellboreModel

    methods
        %-----------------------------------------------------------------%
        function cf = WellborePhaseUpwindFlag(model)

            cf@StateFunction(model);
            cf = cf.dependsOn('massFlux', 'state');

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function flag = evaluateOnDomain(prop, model, state)
            
            % @@TODO: Suppoert multiphase/multicomponent flow
            vt   = value(state.massFlux);
            nph  = model.getNumberOfPhases();
            flag = repmat({vt >= 0}, 1, nph);
            
        end
        %-----------------------------------------------------------------%
        
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