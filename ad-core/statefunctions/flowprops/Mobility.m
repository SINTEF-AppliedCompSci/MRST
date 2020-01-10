classdef Mobility < StateFunction
    properties
    end
    
    methods
        function gp = Mobility(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'RelativePermeability', 'Viscosity'});
        end
        function mob = evaluateOnDomain(prop, model, state)
            [mu, kr] = prop.getEvaluatedDependencies(state, 'Viscosity', 'RelativePermeability');
            mob = cellfun(@(x, y) x./y, kr, mu, 'UniformOutput', false);
            if isfield(model.fluid, 'tranMultR')
                % Pressure dependent mobility multiplier 
                p = model.getProp(state, 'pressure');
                mult = model.fluid.tranMultR(p);
                mob = cellfun(@(x) x.*mult, mob, 'UniformOutput', false);
            end
            % Check for negative values
            mv = cellfun(@(x) min(value(x)), mob);
            if any(mv < 0)
                if model.verbose > 1
                    warning('Negative mobilities detected! Capping to zero.')
                end
                mob = cellfun(@(x) max(x, 0), mob, 'UniformOutput', false);
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
