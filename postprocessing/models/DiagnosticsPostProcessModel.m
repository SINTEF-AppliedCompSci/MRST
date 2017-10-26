classdef DiagnosticsPostProcessModel < ReservoirModel
    properties
        tracerNames
    end

    methods
        function model = DiagnosticsPostProcessModel(G, rock, fluid, varargin)
            require diagnostics
            model = model@ReservoirModel(G, rock, fluid);
            model = merge_options(model, varargin{:});
            model.stepFunctionIsLinear = true;
        end
        function [state, report] = stepFunction(model, state, state0, dt,...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin)
            state.diagnostics = computeTOFandTracer(state0, model.G, model.rock, 'Wells', drivingForces.W);
            report = model.makeStepReport('Converged', true);
        end
    end
end