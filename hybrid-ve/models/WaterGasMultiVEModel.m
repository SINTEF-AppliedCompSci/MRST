classdef WaterGasMultiVEModel < TwoPhaseWaterGasModel
    % Multi-VE version of the gas-water model from the co2lab module
    properties
        
    end
    
    methods
        
        function model = WaterGasMultiVEModel(G, rock, fluid, varargin)
            model = model@TwoPhaseWaterGasModel(G, rock, fluid, nan, nan, varargin{:});
            model = addCoarseOperatorsMultiVE(model);
            model.useCNVConvergence = true;
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, ...
                drivingForces, ...
                varargin)
            [problem, state] = equationsWaterGasMultiVE(model, state0, state , dt , ...
                drivingForces, ...
                varargin{:});
        end
    end
end