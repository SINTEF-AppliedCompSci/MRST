classdef Transmissibility < StateFunction
    % Transmissibility for internal faces. May include an optional
    % pressure-dependent multiplier from a field in the fluid model.
    properties
        
    end
    
    methods
        function pp = Transmissibility(model)
            pp@StateFunction(model);
            if isfield(model.fluid, 'transMult')
                pp = pp.dependsOn('pressure', 'state');
            end
            pp.label = 'T_f';
            assert(isfield(model.operators, 'T'));
            T = value(model.operators.T);
            assert(all(isfinite(T)))
            if any(T < 0)
                warning('Negative transmissibility in %d interfaces', sum(T < 0));
            end
            pp.outputRange = [0, inf];
        end
        
        function T = evaluateOnDomain(sfn, model, state)
            T = model.operators.T;
            if isfield(model.fluid, 'transMult')
                p = model.getProps(state, 'pressure');
                p = model.operators.faceAvg(p);
                T = model.fluid.transMult(p).*T;
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
