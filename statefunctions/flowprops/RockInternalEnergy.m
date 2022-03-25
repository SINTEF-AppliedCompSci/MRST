classdef RockInternalEnergy < StateFunction
%State function for internal energy in the rock

    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = RockInternalEnergy(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('temperature', 'state');
            if isa(model.rock.CpR, 'function_handle')
                % If heat capacity is given as function handle, it's input
                % argumens are pressure and temperature
                gp = gp.dependsOn('pressure', 'state');
            end
            gp.label = 'u_R';
        end
        
        %-----------------------------------------------------------------%
        function uR = evaluateOnDomain(prop, model, state)
            % Get pressure and temperature
            T = model.getProps(state, 'temperature');
            if isa(model.rock.CpR, 'function_handle')
                % Function of pressure and temperature
                p   = model.getProps(state, 'pressure');
                CpR = prop.evaluateFunctionOnDomainWithArguments(model.rock.CpR, p, T); 
            else
                % Constant
                CpR = model.rock.CpR;
            end
            % Compute internal energy
            uR = CpR.*T;
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