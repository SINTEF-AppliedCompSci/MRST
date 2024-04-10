classdef Temperature < StateFunction
   
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = Temperature(model, varargin)
            gp@StateFunction(model, varargin{:});
            switch model.thermalFormulation
                case 'enthalpy'
                    % Enthalpy is the primary variable
                    gp = gp.dependsOn({'PhasePressures', 'Density', 'Enthalpy'});
                case 'temperature'
                    gp = gp.dependsOn({'temperature'}, 'state');
                    % Temperature is the primary variable - we can get it
                    % directly from state
            end
            gp.label = 'T';
        end
        
        %-----------------------------------------------------------------%
        function T = evaluateOnDomain(prop, model, state)
            switch model.thermalFormulation
                case 'enthalpy'
                    error('Not implemented yet')
                case 'temperature'
                    T = state.T;
            end
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