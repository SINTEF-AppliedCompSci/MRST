classdef PhaseMass < StateFunction
%State function for fluid phase mass
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = PhaseMass(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'Density', 'PoreVolume'}, 'PVTPropertyFunctions');
            gp = gp.dependsOn('s', 'state');
            gp.label = 'M_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function mass = evaluateOnDomain(prop, model, state)
            [rho, pv, s] = model.getProps(state, 'Density', 'PoreVolume', 's');
            nph  = model.getNumberOfPhases();
            mass = cell(1, nph);
            for i = 1:nph
                mass{i} = pv.*rho{i}.*s(:, i);
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