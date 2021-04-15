classdef RvMax < StateFunction
    % Maximum Rv (vaporized oil-gas ratio)
    properties
        rvReduction = 0;
    end
    
    methods
        function gp = RvMax(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('PhasePressures');
            gp.label = 'R_v^{max}';
            gp.outputRange = [0, inf];
        end
        
        function rvSat = evaluateOnDomain(prop, model, state)
            p_phase = prop.getEvaluatedDependencies(state, 'PhasePressures');
            pg = p_phase{model.water + model.oil + model.gas};
            rvSat = prop.evaluateFluid(model, 'rvSat', pg);
            if prop.rvReduction > 0 && isfield(state, 'sMax')
                [sOMax, sO] = model.getProps(state, 'somax', 'so');
                sOMax = max(sOMax, sO);
                factor = (sO + 1e-4)./(sOMax + 1e-4);
                rvSat = rvSat.*(factor.^prop.rvReduction);
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
