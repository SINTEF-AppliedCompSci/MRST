classdef EquilibriumConstantModel < EquationOfStateModel
    % Equilibrium constant EOS model for problems where functions of the
    % type K(p, T, z) is sufficient to describe the phase behavior
    properties
    end
    
    methods
        function model = EquilibriumConstantModel(G, fluid, k_values)
            model = model@EquationOfStateModel(G, fluid);
            assert(iscell(k_values));
            assert(all(cellfun(@(x) isa(x, 'function_handle'), k_values)));
            model.equilibriumConstantFunctions = k_values;
        end
        
        function [Z_L, Z_V, f_L, f_V] = getCompressibilityAndFugacity(model, P, T, x, y, z, Z_L, Z_V, varargin)
            R = 8.3144598;
            rhoL = model.PropertyModel.computeMolarDensity(model, P, x, nan, T, true);
            rhoV = model.PropertyModel.computeMolarDensity(model, P, y, nan, T, false);
            assert(all(P > 0), 'Pressures must be positive!')
            Z_L = P./(rhoL.*R.*T);
            Z_V = P./(rhoV.*R.*T);
            if nargout < 3
                return
            end
            K =  model.evaluateEquilibriumConstants(P, T, z);
            if iscell(x)
                [f_L, f_V] = deal(cell(1, numel(x)));
                for i = 1:numel(x)
                    f_V{i} = y{i};
                    f_L{i} = x{i}.*K{i};
                end
            else
                f_V = y;
                f_L = x.*K;
            end
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@PhysicalModel(model, state0, state, dt, drivingForces);
            state = model.setFlag(state);
            state.x = model.computeLiquid(state);
            state.y = model.computeVapor(state);
            [pureLiquid, pureVapor] = model.getFlag(state);
            state.x(pureVapor, :) = state.components(pureVapor, :);
            state.y(pureLiquid, :) = state.components(pureLiquid, :);
            [state.Z_L, state.Z_V] = model.getCompressibilityAndFugacity(state.pressure, state.T, state.x, state.y);
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
