classdef BiotCoupledDilatation < StateFunction
    
    methods
        function gp = BiotCoupledDilatation(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'displacement', 'lambda'}, 'state');
            gp.label = '\nabla\cdot u';
        end
        
        function divu = evaluateOnDomain(prop, model, state)
            [u, lambda] = model.getProps(state, 'displacement', 'lambda');
            divop = model.operators.divop;
            extu = vertcat(u, lambda);
            divu = divop(extu);
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
