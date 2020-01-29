classdef PressureGradient < StateFunction
    properties

    end
    
    methods
        function gp = PressureGradient(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('PhasePressures', 'PVTPropertyFunctions');
            gp = gp.dependsOn('pressure', 'state');
        end
        function dp = evaluateOnDomain(prop, model, state)
            act = model.getActivePhases();
            nph = sum(act);
            Grad = model.operators.Grad;
            
            dp = cell(1, nph);
            if model.FlowPropertyFunctions.CapillaryPressure.pcPresent(model)
                % We have different phase pressures, call gradient once for
                % each phase
                p = model.getProp(state, 'PhasePressures');
                for i = 1:nph
                    dp{i} = Grad(p{i});
                end
            else
                % There is no capillary pressure and a single gradient for
                % the unique pressure is sufficient
                p = model.getProp(state, 'pressure');
                [dp{:}] = deal(Grad(p));
            end
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
