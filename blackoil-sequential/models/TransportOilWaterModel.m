classdef TransportOilWaterModel < TwoPhaseOilWaterModel
    % Two phase oil/water system without dissolution
    properties
        conserveWater
        conserveOil
        staticUpwind
    end
    
    methods
        function model = TransportOilWaterModel(G, rock, fluid, varargin)
            
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            
            model.conserveWater = false;
            model.conserveOil   = true;
            
            model.staticUpwind  = false;

            model = merge_options(model, varargin{:});
            
            assert(~(model.conserveWater && model.conserveOil), ... 
                            'Sequential form only conserves n-1 phases');
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = transportEquationOilWater(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            'solveForOil',   model.conserveOil, ...
                            'solveForWater', model.conserveWater, ...
                            varargin{:});
            
        end
    end
end
