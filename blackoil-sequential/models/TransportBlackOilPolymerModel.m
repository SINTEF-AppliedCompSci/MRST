classdef TransportBlackOilPolymerModel < ThreePhaseBlackOilPolymerModel
    % Two phase oil/water system with polymer
    properties
        conserveWater
        conserveOil
        conserveGas
        upwindType
    end
    
    methods
        function model = TransportBlackOilPolymerModel(G, rock, fluid, varargin)
            
            model = model@ThreePhaseBlackOilPolymerModel(G, rock, fluid);
            
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
        
        function [problem, state] = getEquations(model, state0, state, ...
                dt, drivingForces, varargin)
            [problem, state] = transportEquationBlackOilPolymer(...
                state0, state, model, dt, drivingForces,...
                'solveForWater', model.conserveWater, ...
                'solveForOil'  , model.conserveOil  , ...
                'solveForGas'  , model.conserveGas  , ...
                varargin{:});
        end

    end
end