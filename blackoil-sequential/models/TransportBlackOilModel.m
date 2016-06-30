classdef TransportBlackOilModel < ThreePhaseBlackOilModel
    % Two phase oil/water system without dissolution
    properties
        conserveWater
        conserveOil
        conserveGas
        upwindType
    end
    
    methods
        function model = TransportBlackOilModel(G, rock, fluid, varargin)
            
            model = model@ThreePhaseBlackOilModel(G, rock, fluid);
            
            model.conserveWater = false;
            model.conserveOil   = true;
            model.conserveGas   = true;
            model.upwindType  = 'potential';
            model.useCNVConvergence = false;
            model.nonlinearTolerance = 1e-3;

            model = merge_options(model, varargin{:});
            
            assert(~(model.conserveWater && model.conserveOil && model.conserveGas), ... 
                            'Sequential form only conserves n-1 phases');
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = transportEquationBlackOil(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            'solveForOil',   model.conserveOil, ...
                            'solveForWater', model.conserveWater, ...
                            'solveForGas',   model.conserveGas, ...
                            varargin{:});
            
        end
    end
end
