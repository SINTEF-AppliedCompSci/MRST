classdef ComponentPhaseVelocityFractionalFlowSimpleDG < StateFunction
    properties
    end
    
    methods
        function cf = ComponentPhaseVelocityFractionalFlowSimpleDG(model)
            cf@StateFunction(model);
            cf = cf.dependsOn({'TotalVelocity'});
            assert(isfield(model.fluid, 'f_w'), 'Fractional flow must be specified for water');
        end

        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            
            vT = prop.getEvaluatedDependencies(state, 'TotalVelocity');
            
            assert(ncomp == nph, 'Immiscible assumption');
            
            sw  = model.getProp(state,'sw');
            f_w = model.fluid.f_w(sw);
            
            f_o = 1 - f_w;
            
            v = cell(ncomp, nph);
            v{1, 1} = f_w.*vT;
            v{2, 2} = f_o.*vT;

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
