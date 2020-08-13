classdef RsMax < StateFunction
    % Maximum Rs (dissolved gas-oil ratio)
    properties
    end
    
    methods
        function gp = RsMax(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('pressure', 'state');
            gp.label = 'R_s^{max}';
        end
        
        function rsSat = evaluateOnDomain(prop, model, state)
            p = model.getProp(state, 'pressure');
            if model.disgas
                rsSat = prop.evaluateFluid(model, 'rsSat', p);
                v = value(rsSat);
                if any(v < 0)
                    if model.verbose > 1
                        warning('Negative RsMax detected in %d points', sum(v < 0));
                    end
                    rsSat = max(rsSat, 0);
                end
            else
                rsSat = 0*p;
            end
        end
    end
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
