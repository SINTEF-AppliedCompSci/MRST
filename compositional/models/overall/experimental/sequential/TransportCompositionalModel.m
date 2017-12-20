classdef TransportCompositionalModel < ThreePhaseCompositionalModel
    % Two phase oil/water system without dissolution
    properties
        conserveWater
        staticUpwind
        upwindType
    end
    
    methods
        function model = TransportCompositionalModel(G, rock, fluid, compfluid, varargin)
            
            model = model@ThreePhaseCompositionalModel(G, rock, fluid, compfluid);
            
            model.conserveWater = true;
            model.staticUpwind = false;
            model.upwindType  = 'potential';

            model = merge_options(model, varargin{:});
            if model.water
                model.saturationVarNames = {'sw', 'so', 'sg'};
            else
                model.saturationVarNames = {'so', 'sg'};
            end
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = transportEquationCompositional(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            'solveForWater', model.conserveWater, ...
                            varargin{:});
            
        end
    end
end
