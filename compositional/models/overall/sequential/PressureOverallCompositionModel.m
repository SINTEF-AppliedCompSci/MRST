classdef PressureOverallCompositionModel < OverallCompositionCompositionalModel
    % Two phase oil/water system without dissolution
    properties
        useIncTolPressure
        singlePhaseStrategy = 'numerical';
        twoPhaseStrategy = 'numerical';
        singlePhaseDifferentiation = 'numerical';
        twoPhaseDifferentiation = 'numerical';
    end
    
    methods
        function model = PressureOverallCompositionModel(G, rock, fluid, compFluid, varargin)
            
            model = model@OverallCompositionCompositionalModel(G, rock, fluid, compFluid);
            model.useIncTolPressure = true;
            model = merge_options(model, varargin{:});
            model.outputFluxes = true;
            model.useIncTolComposition = false;
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            state.s = state.s./sum(state.s, 2);
            [problem, state] = equationsCompositional(state0, state, model, dt, ...
                        drivingForces, 'pressure', true, varargin{:});
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            [state, report] = updateState@OverallCompositionCompositionalModel(model, state, problem, dx, drivingForces);
        end

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces) %#ok
            [state, report] = updateAfterConvergence@OverallCompositionCompositionalModel(model, state0, state, dt, drivingForces);
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
