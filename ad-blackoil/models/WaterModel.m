classdef WaterModel < ReservoirModel
    % Two phase oil/water system without dissolution
    properties

    end
    
    methods
        function model = WaterModel(G, rock, fluid, varargin)
            
            model = model@ReservoirModel(G, rock, fluid);
            
            % This is the model parameters for oil/water
            model.oil = false;
            model.gas = false;
            model.water = true;
            %model.addflux=true;
            
            % Blackoil -> use CNV style convergence 
            model.useCNVConvergence = false;
            
            model.saturationVarNames = {'sw'};
            model.wellVarNames = {'qWs', 'bhp'};
            
            model = merge_options(model, varargin{:});
            
            % Setup operators
            %model = model.setupOperators(G, rock, 'deck', model.inputdata);
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsWater(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            % Parent class handles almost everything for us
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
            
            % Update wells based on black oil specific properties
            saturations = model.saturationVarNames;
            wi = strcmpi(saturations, 'sw');
            oi = strcmpi(saturations, 'so');
            gi = strcmpi(saturations, 'sg');

            W = drivingForces.W;
            state.wellSol = assignWellValuesFromControl(model, state.wellSol, W, wi, oi, gi);

        end
    end
end
