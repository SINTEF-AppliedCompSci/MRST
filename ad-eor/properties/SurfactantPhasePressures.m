classdef SurfactantPhasePressures < StateFunction
    properties
    end
    
    methods
        function gp = SurfactantPhasePressures(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('CapillaryPressure', 'FlowPropertyFunctions');
            gp = gp.dependsOn('pressure', 'state');
        end
        
        function p_phase = evaluateOnDomain(prop, model, state)
            p = model.getProp(state, 'Pressure');
            pc = model.getProp(state, 'CapillaryPressure');
            nph = numel(pc);
            p_phase = cell(1, nph);
            for i = 1:nph
                if isempty(pc{i})
                    p_phase{i} = p;
                else
                    p_phase{i} = p + pc{i};
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
