classdef PerforationMobility < StateFunction
    % Mobility in each perforated cell of a well
    properties
    end
    
    methods
        function gp = PerforationMobility(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'FacilityWellMapping'});
            gp = gp.dependsOn({'Mobility'}, 'FlowPropertyFunctions');
            gp.label = '\lambda_{wc}'; 
        end
        
        function mobw = evaluateOnDomain(prop, model, state)
            rm = model.ReservoirModel;
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            mob = rm.FlowPropertyFunctions.get(rm, state, 'Mobility');
            if ~iscell(mob)
                mob = {mob};
            end
            mobw = applyFunction(@(x) x(map.cells), mob);
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
