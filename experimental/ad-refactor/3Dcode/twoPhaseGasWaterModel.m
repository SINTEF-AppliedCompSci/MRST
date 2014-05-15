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
        
        function state = updateState(model, state, problem, dx, drivingForces) %#ok
            
            dp      = dx{problem.indexOfPrimaryVarible('pressure')};
            ds      = dx{problem.indexOfPrimaryVarible('sG')};
            dqWs    = dx{problem.indexOfPrimaryVarible('qWs')};
            dqGs    = dx{problem.indexOfPrimaryVarible('qGs')};
            dpBHP   = dx{problem.indexOfPrimaryVarible('bhp')};
            
            dsMax   = model.dsMax; 
            dpMax   = model.dpMax .* state.pressure;
            dMaxBHP = model.dpMax .* vertcat(state.wellSol.bhp);

            ds      = truncate(-dsMax   , ds                , dsMax   );
            dp      = truncate(-dpMax   , dp                , dpMax   );
            dpBHP   = truncate(-dMaxBHP , dpBHP             , dMaxBHP );
            sG      = truncate(0        , state.s(:,2) + ds , 1       );
            
            state.pressure = state.pressure + dp;
            state.s = [1-sG, sG];

            for w = 1:numel(state.wellSol)
                state.wellSol(w).bhp = state.wellSol(w).bhp + dpBHP(w);
                state.wellSol(w).qWs = state.wellSol(w).qWs + dqWs(w);
                state.wellSol(w).wGs = state.wellSol(w).qGs + dqGs(w);
            end
        end
    end
    
end
% ============================================================================

function res = truncate(m, val, M)
% ensure all components in 'val' are in the range [m, M]
    res = max(m, min(val, M)); 
end

function t = computeTemperatureField(G, tsurf, tgrad)
    t = tsurf * ones(G.cells.num, 1) + G.cells.centroids(:,3) * tgrad / 1000;
end


