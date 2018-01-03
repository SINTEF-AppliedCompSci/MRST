classdef TransportNaturalVariablesModel < NaturalVariablesCompositionalModel
    properties
        staticUpwind
        upwindType
    end
    
    methods
        function model = TransportNaturalVariablesModel(G, rock, fluid, compFluid, varargin)
            
            model = model@NaturalVariablesCompositionalModel(G, rock, fluid, compFluid);
            model.staticUpwind = false;
            model.upwindType  = 'potential';
            model = merge_options(model, varargin{:});
            model.EOSModel.fastDerivatives = false;
            model.allowLargeSaturations = true;
            model.maxPhaseChangesNonLinear = inf;
            if model.water
                model.saturationVarNames = {'sw', 'so', 'sg'};
            else
                model.saturationVarNames = {'so', 'sg'};
            end
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = transportEquationsNaturalVars(state0, state, model, dt, ...
                            drivingForces, varargin{:});
        end
    end
end


