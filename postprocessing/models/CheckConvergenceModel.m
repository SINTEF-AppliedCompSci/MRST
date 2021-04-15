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
            [problem, state] = model.parentModel.getEquations(state0, state, dt, drivingForces, 'resOnly', true, 'iteration', inf);
            [convergence, values, names] = model.parentModel.checkConvergence(problem);
            report = model.makeStepReport(...
                                    'Failure',         false, ...
                                    'Converged',       true, ...
                                    'ResidualsConverged', true(size(convergence)), ...
                                    'Residuals',       values ...
                                    );
            state = model.parentModel.reduceState(state, false);
            state = model.parentModel.updateAfterConvergence(state0, state, dt, drivingForces);
            
            state.convergenceStatus = struct('Residuals', values, 'ResidualsConverged', convergence,...
                                             'Converged', all(convergence), 'Names', {names'});
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

        function [model, state] = prepareReportstep(model, state, varargin)
            model.parentModel = model.parentModel.prepareReportstep(state, varargin{:});
        end

        function [model, state] = prepareTimestep(model, state, varargin)
            model.parentModel = model.parentModel.prepareReportstep(state, varargin{:});
        end
    end
    
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
