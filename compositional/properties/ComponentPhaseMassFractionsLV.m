classdef ComponentPhaseMassFractionsLV < StateFunction
    % Component mass fractions
    properties
    end
    
    methods
        function gp = ComponentPhaseMassFractionsLV(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'ComponentPhaseMoleFractions'});
            gp.label = 'X_{i,\alpha}';
        end

        function mass = evaluateOnDomain(prop, model, state)
            eos = model.EOSModel;
            moles = prop.getEvaluatedDependencies(state, 'ComponentPhaseMoleFractions');
            
            mass = moles;
            isEoS = model.getEoSComponentMask();
            for i = 1:size(mass, 2)
                mass(isEoS, i) = eos.getMassFraction(moles(isEoS, i));
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
