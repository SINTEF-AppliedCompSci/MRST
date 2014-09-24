classdef PressureOilWaterModel < TwoPhaseOilWaterModel
    % Two phase oil/water system without dissolution
    properties

    end
    
    methods
        function model = PressureOilWaterModel(G, rock, fluid, varargin)
            
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);

            model = merge_options(model, varargin{:});

            % Ensure simple tolerances
            model.useCNVConvergence = false;
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = pressureEquationOilWater(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            
        end
    end
end
