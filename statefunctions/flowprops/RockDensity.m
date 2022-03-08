classdef RockDensity < StateFunction
%State function for rock density
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = RockDensity(model, varargin)
            gp@StateFunction(model, varargin{:});
            if isa(model.rock.rhoR, 'function_handle')
                % If density is given as function handle, it's input
                % argumens are pressure and temperature
                gp = gp.dependsOn({'pressure', 'temperature'}, 'state');
            end
            gp.label = '\rho_R';
        end
        
        %-----------------------------------------------------------------%
        function rhoR = evaluateOnDomain(prop, model, state)
            if isa(model.rock.rhoR, 'function_handle')
                % Function of pressure and temperature
                [p, T] = model.getProps(state, 'pressure', 'temperature');
                rhoR = model.rock.rhoR(p,T);
            else
                % Constant
                rhoR = model.rock.rhoR;
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