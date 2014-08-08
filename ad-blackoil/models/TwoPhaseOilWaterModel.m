classdef TwoPhaseOilWaterModel < ReservoirModel
    % Two phase oil/water system without dissolution
    properties

    end
    
    methods
        function model = TwoPhaseOilWaterModel(G, rock, fluid, varargin)
            
            model = model@ReservoirModel(G, rock, fluid);
            
            % This is the model parameters for oil/water
            model.oil = true;
            model.gas = false;
            model.water = true;
            
            model.saturationVarNames = {'sw', 'so'};
            model.wellVarNames = {'qWs', 'qOs', 'bhp'};
            
            model = merge_options(model, varargin{:});
            
            % Setup operators
            model = model.setupOperators(G, rock, 'deck', model.inputdata);
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsOilWater(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});
            
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            state = updateStateBlackOilGeneric(model, state, problem, dx, drivingForces);
            report = struct();
        end
    end
end
