classdef WellIndex < StateFunction
    % Get the well index (productivity or injectivity index, equivialent to
    % a transmissibility between the reservoir and the wellbore).
    properties

    end
    
    methods
        function gp = WellIndex(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('FacilityWellMapping');
            gp.label = 'W_I';
        end
        function WI = evaluateOnDomain(prop, model, state)
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            wellSol = state.wellSol;
            W = model.getWellStruct(map.active);
            if isempty(W)
                WI = [];
            else
                cstatus = vertcat(wellSol(map.active).cstatus);
                WI = vertcat(W.WI).*cstatus;
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
