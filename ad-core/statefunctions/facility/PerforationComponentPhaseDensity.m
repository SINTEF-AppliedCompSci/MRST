classdef PerforationComponentPhaseDensity < StateFunction
    % Component density to used for each well connection
    properties

    end
    
    methods
        function gp = PerforationComponentPhaseDensity(varargin)
            gp = gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'ComponentPhaseDensity'}, 'FlowPropertyFunctions');
            gp = gp.dependsOn({'FacilityWellMapping'});
            gp.label = '\rho_{wc}';
        end
        
        function rhoc = evaluateOnDomain(prop, model, state)
            rm = model.ReservoirModel;
            fpr = rm.FlowPropertyFunctions;
            rhoc = fpr.get(rm, state, 'ComponentPhaseDensity');
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            for i = 1:numel(rhoc)
                if ~isempty(rhoc{i})
                    rhoc{i} = rhoc{i}(map.cells);
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
