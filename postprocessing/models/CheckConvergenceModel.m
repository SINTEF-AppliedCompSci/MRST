classdef CheckConvergenceModel < PhysicalModel
    properties
        parentModel
    end
    methods
        function model = CheckConvergenceModel(parent)
            assert(isa(parent, 'PhysicalModel'));
            model = model@PhysicalModel(parent.G);
            model.parentModel = parent;
        end
        function [state, report] = stepFunction(model, state, state0, dt,...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin)
            problem = model.parentModel.getEquations(state0, state, dt, drivingForces, 'resOnly', true, 'iteration', inf);
            [convergence, values, names] = model.parentModel.checkConvergence(problem);
            report = model.makeStepReport(...
                                    'Failure',         false, ...
                                    'Converged',       true, ...
                                    'ResidualsConverged', true(size(convergence)), ...
                                    'Residuals',       values ...
                                    );
            state.convergenceStatus = struct('Residuals', values, 'ResidualsConverged', convergence,...
                                             'Converged', all(convergence), 'Names', {names});
        end
        
        function varargout = getActivePhases(model)
            varargout = cell(1, nargout);
            [varargout{:}] = model.parentModel.getActivePhases();
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
            return
        end

        function [fn, index] = getVariableField(model, name)
            [fn, index] = model.parentModel.getVariableField(name);
        end

    end
    
end