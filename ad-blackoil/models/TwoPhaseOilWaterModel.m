classdef TwoPhaseOilWaterModel < ThreePhaseBlackOilModel
% Two phase oil/water system without dissolution
properties

end

methods
    function model = TwoPhaseOilWaterModel(G, rock, fluid, varargin)
        model = model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});

        % This is the model parameters for oil/water
        model.oil = true;
        model.gas = false;
        model.water = true;

        % Blackoil -> use CNV style convergence 
        model.useCNVConvergence = true;

        model.saturationVarNames = {'sw', 'so'};
        model.wellVarNames = {'qWs', 'qOs', 'bhp'};

        model = merge_options(model, varargin{:});
    end
    
    % --------------------------------------------------------------------%
    function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        [problem, state] = equationsOilWater(state0, state, model,...
                        dt, ...
                        drivingForces,...
                        varargin{:});

    end
end
end
