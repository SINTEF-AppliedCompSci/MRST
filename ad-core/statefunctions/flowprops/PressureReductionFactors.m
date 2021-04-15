classdef PressureReductionFactors < StateFunction
    % Component weighting factors used to form a pressure equation -
    % immiscible base class
    properties
    end
    
    methods
        function gp = PressureReductionFactors(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'Density', 'PoreVolume'});
            gp.label = 'w_i^p';
        end

        function w = evaluateOnDomain(prop, model, state)
            [rho, pv] = prop.getEvaluatedDependencies(state, 'Density', 'PoreVolume');
            nph = numel(rho);
            w = cell(1, nph);
            for ph = 1:nph
                w{ph} = pv./rho{ph}; % Scale final pressure equation with pore-volume
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
