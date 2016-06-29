classdef PressureOilWaterModel < TwoPhaseOilWaterModel
    % Two phase oil/water system without dissolution
    properties
        incTolPressure
        useIncTol
    end
    
    methods
        function model = PressureOilWaterModel(G, rock, fluid, varargin)
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            
            model.incTolPressure = 1e-3;
            model.useIncTol = true;
            model = merge_options(model, varargin{:});

            % Ensure simple tolerances
            model.useCNVConvergence = false;
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = pressureEquationOilWater(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            p0 = state.pressure;
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
            range = max(p0) - min(p0);
            if range == 0
                range = 1;
            end
            state.dpRel = (state.pressure - p0)./range;
        end
        
        function  [convergence, values, names] = checkConvergence(model, problem)
            [convergence, values, names] = checkConvergence@PhysicalModel(model, problem);
            if ~isnan(problem.iterationNo) && model.useIncTol
                if problem.iterationNo  > 1
                    values(1) = norm(problem.state.dpRel, inf);
                    convergence = all([values(1) < model.incTolPressure, values(2:end) < model.nonlinearTolerance]);
                else
                    values(1) = inf;
                    convergence = false;
                end
                names{1} = 'Delta P';
            end
        end
    end
end
