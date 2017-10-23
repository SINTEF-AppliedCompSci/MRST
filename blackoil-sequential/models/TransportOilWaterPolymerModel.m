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
            model.useCNVConvergence = false;
            model.nonlinearTolerance = 1e-3;
            
            model = merge_options(model, varargin{:});
            
            assert(~(model.conserveWater && model.conserveOil), ... 
                            'Sequential form only conserves n-1 phases');

        end
        
        function [problem, state] = getEquations(model, state0, state, ...
                dt, drivingForces, varargin)
            [problem, state] = transportEquationOilWaterPolymer(...
                state0, state, model, dt, drivingForces,...
                'solveForOil',   model.conserveOil, ...
                'solveForWater', model.conserveWater, ...
                varargin{:});
        end
        
    end
end
