classdef TransportBlackOilModel < ThreePhaseBlackOilModel
    % Two phase oil/water system without dissolution
    properties
        conserveWater
        conserveOil
        conserveGas
        staticUpwind
    end
    
    methods
        function model = TransportBlackOilModel(G, rock, fluid, varargin)
            
            model = model@ThreePhaseBlackOilModel(G, rock, fluid);
            
            model.conserveWater = false;
            model.conserveOil   = true;
            model.conserveGas   = true;
            
            model.staticUpwind  = false;

            model = merge_options(model, varargin{:});
            
            assert(~(model.conserveWater && model.conserveOil && model.conserveGas), ... 
                            'Sequential form only conserves n-1 phases');
            
            % Ensure simple tolerances
            model.useCNVConvergence = false;
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
        function [convergence, values] = checkConvergence(model, problem, varargin)
            [convergence, values] = checkConvergence@PhysicalModel(model, problem, varargin{:});
            % Always make at least one update so that the problem actually changes.
            convergence = convergence && problem.iterationNo > 1;
        end
    end
end
