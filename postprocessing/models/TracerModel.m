classdef TracerModel < ReservoirModel
    properties
        tracerNames
    end

    methods
        function model = TracerModel(G, rock, fluid, varargin)
            model = model@ReservoirModel(G, rock, fluid);
            model = merge_options(model, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsTracer(state0, state, model, dt, ...
                            drivingForces, varargin{:});

        end

        function n = getNumberOfTracers(model)
            n = numel(model.tracerNames);
        end

        function state = validateState(model, state)
            % Check parent class
            state = validateState@ReservoirModel(model, state);
            nc = model.G.cells.num;
            if isfield(state, 'tracer')
                model.checkProperty(state, 'tracer', nc, 1);
            else
                nt = model.getNumberOfTracers();
                state.tracer = zeros(nc, nt);
            end
        end

        function [fn, index] = getVariableField(model, name)
            tsub = strcmpi(model.tracerNames, name);
            if any(tsub)
                fn = 'tracer';
                index = find(tsub);
            else
                switch(lower(name))
                    case 'tracer'
                        fn = 'tracer';
                        index = ':';
                    otherwise
                        % Basic phases are known to the base class
                        [fn, index] = getVariableField@ReservoirModel(model, name);
                end
            end
        end
    end
end