classdef PressureOverallCompositionModel < OverallCompositionCompositionalModel
    % Two phase oil/water system without dissolution
    properties
        useIncTolPressure
    end
    
    methods
        function model = PressureOverallCompositionModel(G, rock, fluid, compFluid, varargin)
            
            model = model@OverallCompositionCompositionalModel(G, rock, fluid, compFluid);
            model.useIncTolPressure = true;            
            model = merge_options(model, varargin{:});
            model.outputFluxes = true;
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            state.s = state.s./sum(state.s, 2);
            [problem, state] = equationsCompositional(state0, state, model, dt, ...
                        drivingForces, 'pressure', true, varargin{:});
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            [state, report] = updateState@OverallCompositionCompositionalModel(model, state, problem, dx, drivingForces);
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces) %#ok
            [state, report] = updateAfterConvergence@OverallCompositionCompositionalModel(model, state0, state, dt, drivingForces);
        end

    end
end
