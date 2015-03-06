classdef TransportOilWaterPolymerModel < OilWaterPolymerModel
    % Two phase oil/water system with polymer
    properties
        conserveWater
        conserveOil
        staticUpwind
    end
    
    methods
        function model = TransportOilWaterPolymerModel(G, rock, fluid, varargin)
            
            model = model@OilWaterPolymerModel(G, rock, fluid);
            
            model.conserveWater = true;
            model.conserveOil   = false;
            
            model.staticUpwind  = false;

            model = merge_options(model, varargin{:});
            
            assert(~(model.conserveWater && model.conserveOil), ... 
                            'Sequential form only conserves n-1 phases');
            
            % Ensure simple tolerances
            model.useCNVConvergence = false;
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = transportEquationOilWaterPolymer(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            'solveForOil',   model.conserveOil, ...
                            'solveForWater', model.conserveWater, ...
                            varargin{:});
            
        end
        
        function [convergence, values] = checkConvergence(model, problem, varargin)
            [convergence, values] = checkConvergence@PhysicalModel(model, problem, varargin{:});
            % Always make at least one update so that the problem actually changes.
            convergence = convergence && problem.iterationNo > 1;
        end
    end
end
