classdef MpfaKgrad < StateFunction
    
    methods
        function pp = MpfaKgrad(model)
            pp@StateFunction(model);
            pp = pp.dependsOn('PhasePressures', 'PVTPropertyFunctions');
        end
        
        function v = evaluateOnDomain(prop, model, state)
            p_phase = model.getProp(state, 'PhasePressures');
            fluxop = model.operators.fluxop;
            nph = numel(p_phase);
            v = cell(1, nph);
            for i = 1 : nph
                v{i} = fluxop(p_phase{i});
            end
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

