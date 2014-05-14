classdef twoPhaseGasWaterModel < physicalModel

    properties
        % (constant) temperature field
        t
        
    end
    
    % ============================================================================
    methods
        % ------------------------------------------------------------------------
        function model = twoPhaseGasWaterModel(G, rock, fluid, tsurf, tgrad, varargin)
            
            model.G     = G;
            model.fluid = fluid;
            model.oil   = false;
            model.gas   = true;
            model.water = true;
            model.t     = computeTemperatureField(G, tsurf, tgrad);
            model.name  = 'GasWater_2ph';
            model       = model.setupOperators(G, rock);
        end
        % ------------------------------------------------------------------------
        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            [problem, state] = equationsGasWater(state0, state , dt , ...
                                                 model.G         ,    ...
                                                 drivingForces   ,    ...
                                                 model.operators ,    ...
                                                 model.fluid     ,    ...
                                                 model.t         ,    ...
                                                 varargin{:});
        end
        % ------------------------------------------------------------------------
        function state = updateState(model, state, problem, dx, drivingForces)
            error('unimplemneted');
        end
    end
    
end
% ============================================================================

function t = computeTemperatureField(G, tsurf, tgrad)
    t = tsurf * ones(G.cells.num, 1) + G.cells.centroids(:,3) * tgrad / 1000;
end


