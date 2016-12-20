classdef TransportOilWaterPolymerModel < OilWaterPolymerModel
    % Two phase oil/water system with polymer
    properties
        conserveWater
        conserveOil
        upwindType
    end
    
    methods
        function model = TransportOilWaterPolymerModel(G, rock, fluid, varargin)
            
            model = model@OilWaterPolymerModel(G, rock, fluid);
            
            model.conserveWater = false;
            model.conserveOil   = true;
            
            model.upwindType  = 'potential';

            model = merge_options(model, varargin{:});
            
            assert(~(model.conserveWater && model.conserveOil), ... 
                            'Sequential form only conserves n-1 phases');
            
            % Ensure simple tolerances
            model.useCNVConvergence = false;
        end
        
        function [problem, state] = getEquations(model, state0, state, ...
                dt, drivingForces, varargin)
            [problem, state] = transportEquationOilWaterPolymer(...
                state0, state, model, dt, drivingForces,...
                'solveForOil',   model.conserveOil, ...
                'solveForWater', model.conserveWater, ...
                varargin{:});
        end
        
        function [convergence, values, names] = checkConvergence(model, ...
                problem, varargin)
            
            % TODO: HACK
            % To use the convergence methods, we need to temporarily remove
            % the phase from the model that we do not have an equation for.
            if model.conserveWater
                model.oil = false;
            elseif model.conserveOil
                model.water = false;
            end
            [convergence, values, names] = ...
                checkConvergence@OilWaterPolymerModel(model, problem, ...
                varargin{:});
            if model.conserveWater
                model.oil = true;
            elseif model.conserveOil
                model.water = true;
            end
            
            % Always make at least one update so that the problem actually 
            % changes.
            convergence = convergence && problem.iterationNo > 1;
        end
    end
end
