classdef Strain < StateFunction
    
    methods
        function gp = Strain(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('displacement', 'state');
            gp.label = 'strain';
        end
        
        function strain = evaluateOnDomain(prop, model, state)
            op = model.operators.global_strain;
            u = model.getProp(state, 'displacement');
            error('not yet implemented for mpsa')
            strain = op*u;
        end
    end
end

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}

