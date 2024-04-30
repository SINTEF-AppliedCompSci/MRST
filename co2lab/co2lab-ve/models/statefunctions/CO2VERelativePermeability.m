classdef CO2VERelativePermeability < BaseRelativePermeability

    properties
    end
    
    methods
        function gp = CO2VERelativePermeability(model, varargin)
            gp@BaseRelativePermeability(model, varargin{:});
            
            if model.hysteresis
                gp = gp.dependsOn({'sGmax', 'pressure'}, 'state');
            end
        end

        function kr = evaluateOnDomain(prop, model, state)
            sG = model.getProp(state, 'sg');
            
            if model.hysteresis
                sGmax = model.getProp(state, 'sGmax');
            else
                sGmax = sG;
            end
            sW = 1 - sG;
            p = model.getProp(state, 'pressure');
            krW = prop.evaluateFluid(model, 'krW', sW, p, 'sGmax', sGmax);
            krG = prop.evaluateFluid(model, 'krG', sG, p, 'sGmax', sGmax);
            kr = {krW, krG};
        end
    end
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
