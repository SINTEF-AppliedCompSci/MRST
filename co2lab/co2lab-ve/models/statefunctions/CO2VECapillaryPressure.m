classdef CO2VECapillaryPressure < StateFunction
   
    properties
    end
    
    methods
        function prop = CO2VECapillaryPressure(model, varargin)
            prop = prop@StateFunction(model, varargin{:});
            prop = prop.dependsOn({'s', 'sGmax', 'pressure'}, 'state');
            prop.label = 'p_{c}';
        end
        
        function pc = evaluateOnDomain(prop, model, state)
            sG = model.getProp(state, 'sg');
            sGmax = model.getProp(state, 'sGmax');
            pressure = model.getProp(state, 'pressure');
            pc = prop.evaluateFluid(model, 'pcWG', sG, pressure, 'sGmax', sGmax);
            pc = {0*pc, pc}; % we need one column for water, another for CO2,
                             % indicating what to be added to the reference 
                             % pressure to get the phase pressure.  Since
                             % reference pressure is the same as water
                             % pressure here, the corresponding column is zero.
        end
        
        function res = pcPresent(prop, model)
            res = true;
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
