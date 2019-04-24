classdef WrapperModel < PhysicalModel
    % Wrapper model which can be inherited for operations on the parent
    % model.
    properties
        parentModel
    end
    
    methods
        function model = WrapperModel(parentModel, varargin)
            model = model@PhysicalModel(parentModel.G);
            model = merge_options(model, varargin{:});
            
            model.parentModel = parentModel;
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)
            [vars, names, origin] = model.parentModel.getPrimaryVariables(state);
        end

        function state = initStateAD(model, state, vars, names, origin)
            state = model.parentModel.initStateAD(state, vars, names, origin);
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = model.parentModel.getEquations(state0, state, dt, drivingForces, varargin{:});
        end
        
        function forces = getValidDrivingForces(model)
            forces = model.parentModel.getValidDrivingForces();
        end

        function state = validateState(model, state)
            state = model.parentModel.validateState(state);
        end

        function [model, state] = updateForChangedControls(model, state, forces)
            [model.parentModel, state] = model.parentModel.updateForChangedControls(state, forces);
        end

        function model = validateModel(model, varargin)
            model.parentModel = model.parentModel.validateModel(varargin{:});
        end

        function [fn, index] = getVariableField(model, name, varargin)
            [fn, index] = model.parentModel.getVariableField(name, varargin{:});
        end

        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            % Prepare state and model (temporarily) before solving a time-step
            [model.parentModel, state] = model.parentModel.prepareTimestep(state, state0, dt, drivingForces);
        end

        function [model, state] = prepareReportstep(model, state, state0, dt, drivingForces)
            % Prepare state and model (temporarily) before solving a report-step
            [model.parentModel, state] = model.parentModel.prepareReportstep(state, state0, dt, drivingForces);
        end
    end
end