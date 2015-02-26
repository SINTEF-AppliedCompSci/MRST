classdef PressureOilWaterPolymerModel < OilWaterPolymerModel
    % Two phase oil/water system with polymer
    properties

    end
    
    methods
        function model = PressureOilWaterPolymerModel(G, rock, fluid, varargin)
            
            model = model@OilWaterPolymerModel(G, rock, fluid);

            model = merge_options(model, varargin{:});

            % Ensure simple tolerances
            model.useCNVConvergence = false;
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = pressureEquationOilWaterPolymer(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            
        end
        function [convergence, values] = checkConvergence(model, problem, varargin)
            [convergence, values] = checkConvergence@PhysicalModel(model, problem, varargin{:});
            % Always make at least one update so that the problem actually changes.
            convergence = convergence && problem.iterationNo > 1;
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
            % state = updateWellCellSaturationsExplicit(model, state, problem, dx, drivingForces);
        end
    end
end
