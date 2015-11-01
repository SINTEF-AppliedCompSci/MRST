classdef twoPhaseGasWaterModel < ReservoirModel

    properties
        % (constant) temperature field
        t
        name
    end
    
    % ============================================================================
    methods
        % ------------------------------------------------------------------------
        function model = twoPhaseGasWaterModel(G, rock, fluid, tsurf, tgrad, varargin)
           
            model = model@ReservoirModel(G); 
           
            %model.G     = G;
            model.fluid = fluid;
            model.oil   = false;
            model.gas   = true;
            model.water = true;
            model.t     = computeTemperatureField(G, tsurf, tgrad);
            model.name  = 'GasWater_2ph';
            model.saturationVarNames = {'sw', 'sg'}; % @@ Design: ideally, we
                                                     % should not have to set
                                                     % both this _and_ the
                                                     % 'water', 'gas' and 'oil'
                                               % flags above. Check with
                                               % maintainer of parent class.
            model.wellVarNames = {'qWs', 'qGs', 'bhp'};
            model.gravity = gravity;
            model       = model.setupOperators(G, rock);
        end

    end
    
    methods %(Access = protected)

        % ------------------------------------------------------------------------
        
        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            [problem, state] = equationsGasWater(model, state0, state , dt , ...
                                                 drivingForces             , ...
                                                 varargin{:});
        end

        % ------------------------------------------------------------------------
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
            
           [state, report] = updateState@ReservoirModel(model, state, problem, dx, ...
                                                        drivingForces);
           state.sGmax = max(state.sGmax,state.s(:,2));
           state.sGmax = min(1,state.sGmax);
           state.sGmax = max(0,state.sGmax);
        end
    end
    
end
% ============================================================================

function res = truncate(m, val, M)
% ensure all components in 'val' are in the range [m, M]
    if ~isempty(val)
        res = max(m, min(val, M)); 
    else
        res = [];
    end
end

function t = computeTemperatureField(G, tsurf, tgrad)
    t = tsurf * ones(G.cells.num, 1) + G.cells.centroids(:,3) * tgrad / 1000;
end


