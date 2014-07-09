classdef TwoPhaseOilWaterModel < PhysicalModel
    % two phase oil / water system
    properties

    end
    
    methods
        function model = TwoPhaseOilWaterModel(G, rock, fluid, varargin)
            
            model = model@PhysicalModel(G, rock, fluid);
            
            % This is the model parameters for oil/water
            model.oil = true;
            model.gas = false;
            model.water = true;
            
            model.componentNames = {'sw', 'so'};
            
            model = merge_options(model, varargin{:});
            
            % Name and operators
            model.name = 'OilWater_2ph';
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
