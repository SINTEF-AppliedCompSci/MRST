classdef PressureCompositionalModel < OverallCompositionCompositionalModel
    % Two phase oil/water system without dissolution
    properties
        useIncTolPressure
        useMassFlux
    end
    
    methods
        function model = PressureCompositionalModel(G, rock, fluid, compFluid, varargin)
            
            model = model@OverallCompositionCompositionalModel(G, rock, fluid, compFluid);
            model.incTolPressure = 1e-3;
            model.useIncTolPressure = true;
            model.dpMaxRel = 0.2;
            model.useMassFlux = true;
            
            model = merge_options(model, varargin{:});
            if model.water
                model.saturationVarNames = {'sw', 'so', 'sg'};
            else
                model.saturationVarNames = {'so', 'sg'};
            end
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            if model.useMassFlux
                [problem, state] = pressureEquationCompositional(state0, state, model, dt, ...
                            drivingForces, varargin{:});
            else
                [problem, state] = equationsCompositional(state0, state, model, dt, ...
                            drivingForces, 'pressure', true, varargin{:});
                state.massFlux = false;
            end

        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            [state, report] = updateState@ThreePhaseCompositionalModel(model, state, problem, dx, drivingForces);
        end
        
        function  [convergence, values, names] = checkConvergence(model, problem)
            [convergence, values, names] = checkConvergence@ThreePhaseCompositionalModel(model, problem);
            if ~isnan(problem.iterationNo) && model.useIncTol
                if problem.iterationNo  > 1
                    values(1) = norm(problem.state.dpRel, inf);
                else
                    values(1) = inf;
                    convergence(1) = false;
                end
                names{1} = 'Delta P';
            end
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces) %#ok
            [state, report] = updateAfterConvergence@ThreePhaseCompositionalModel(model, state0, state, dt, drivingForces);
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
