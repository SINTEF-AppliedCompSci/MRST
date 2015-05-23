classdef PressureOilWaterPolymerModel < OilWaterPolymerModel
    % Two phase oil/water system with polymer
    properties
        
    end
    
    methods
        function model = PressureOilWaterPolymerModel(G, rock, fluid, ...
                varargin)
            
            model = model@OilWaterPolymerModel(G, rock, fluid);

            model = merge_options(model, varargin{:});

            % Ensure simple tolerances
            model.useCNVConvergence = false;
            
        end
        
        function [problem, state] = getEquations(model, state0, state, ...
                dt, drivingForces, varargin)
            [problem, state] = pressureEquationOilWaterPolymer(state0, ...
                state, model, dt, drivingForces, varargin{:});
        end
        
        function [convergence, values, names] = checkConvergence(model, ...
                problem, varargin)
            [convergence, values, names] = ...
                checkConvergence@OilWaterPolymerModel(model, problem, ...
                varargin{:});
            
            % Always make at least one update so that the problem actually 
            % changes.
            convergence = convergence && problem.iterationNo > 1;
        end
        
        function [state, report] = updateAfterConvergence(model, ...
                state0, state, dt, drivingForces)
            [state, report] = ...
                updateAfterConvergence@OilWaterPolymerModel(model, ...
                state0, state, dt, drivingForces);
            if model.polymer
                % Special hack for the sequential solver with shear
                % thinning. See equations for details.
                for w=1:numel(state.wellSol)
                    state.wellSol(w).poly_prev = ...
                        drivingForces.Wells(w).poly;
                end
            end
        end
        
    end
end
