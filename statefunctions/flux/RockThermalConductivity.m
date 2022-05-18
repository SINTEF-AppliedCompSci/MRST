classdef RockThermalConductivity < StateFunction
%State function for computing the thermal conductivity of the rock
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = RockThermalConductivity(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PoreVolume'}, 'PVTPropertyFunctions');
            gp = gp.dependsOn({'pressure', 'temperature'}, 'state');
            gp.label = '\Lambda_R';
        end
        
        %-----------------------------------------------------------------%
        function lambdaR = evaluateOnDomain(prop, model, state)
            [p, T, pv] = model.getProps(state, 'pressure', 'Temperature', 'PoreVolume');
            if model.dynamicHeatTransRock
                lambdaR = prop.evaluateFunctionOnDomainWithArguments(model.rock.lambdaR, p, T);
            else
                lambdaR = model.rock.lambdaR;
            end
            vol     = model.operators.vol;
            lambdaR = lambdaR.*(vol - pv)./vol;
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