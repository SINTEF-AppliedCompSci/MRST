classdef PostProcessModel < PhysicalModel
    properties
        parentModel
        postProcessModels
        postProcessNonLinearSolvers
    end
    
    methods
        function model = PostProcessModel(parentModel, postProcessModels, varargin)
            model = model@PhysicalModel(parentModel.G);
            model = merge_options(model, varargin{:});
            
            model.parentModel = parentModel;
            assert(iscell(postProcessModels));
            model.postProcessModels = postProcessModels;
            if isempty(model.postProcessNonLinearSolvers)
                nls = NonLinearSolver();
                model.postProcessNonLinearSolvers = cell(numel(model.postProcessModels), 1);
                [model.postProcessNonLinearSolvers{:}] = deal(nls);
            end
        end
        
        function [state, report] = stepFunction(model, state, state0, dt,...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin)
            [state, report] = model.parentModel.stepFunction(...
                state, state0, dt, drivingForces, linsolve, ...
                nonlinsolve, iteration, varargin);
            if report.Converged
                % We are done, apply post-processing
                for i = 1:numel(model.postProcessModels)                    
                    submodel = model.postProcessModels{i};
                    nls = model.postProcessNonLinearSolvers{i};
                    
                    arg = submodel.getDrivingForces(drivingForces);
                    state = nls.solveTimestep(state0, dt, submodel, varargin{:}, 'initialGuess', state, arg{:});
                end
            end
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
            for i = 1:numel(model.postProcessModels)
                state = model.postProcessModels{i}.validateState(state);
            end
        end

        function [model, state] = updateForChangedControls(model, state, forces)
            [model.parentModel, state] = model.parentModel.updateForChangedControls(state, forces);
            for i = 1:numel(model.postProcessModels)
                [model.postProcessModels{i}, state] = model.postProcessModels{i}.updateForChangedControls(state, forces);
            end
        end

        function model = validateModel(model, varargin)
            model.parentModel = model.parentModel.validateModel(varargin{:});
            for i = 1:numel(model.postProcessModels)
                model.postProcessModels{i} = model.postProcessModels{i}.validateModel(varargin{:});
            end
        end

        function [fn, index] = getVariableField(model, name)
            [fn, index] = model.parentModel.getVariableField(name);
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
