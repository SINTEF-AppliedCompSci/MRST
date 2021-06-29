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
        
        function schedule = validateSchedule(model, schedule)
            schedule = model.parentModel.validateSchedule(schedule);
        end

        function [state, report] = updateState(model, varargin)
            [state, report] = model.parentModel.updateState(varargin{:});
        end
        
        function [convergence, values, names] = checkConvergence(model, varargin)
            [convergence, values, names] = model.parentModel.checkConvergence(varargin{:});
        end
        
        function [model, state] = updateForChangedControls(model, state, forces)
            [model.parentModel, state] = model.parentModel.updateForChangedControls(state, forces);
        end

        function model = validateModel(model, varargin)
            setDefaults = isempty(model.parentModel.getStateFunctionGroupings());
            model.parentModel = model.parentModel.validateModel(varargin{:});
            model = model.setupStateFunctionGroupings(setDefaults);
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
        
        function [vararg, control] = getDrivingForces(model, control)
            [vararg, control] = model.parentModel.getDrivingForces(control);
        end
        
        function [state, report] = updateAfterConvergence(model, varargin)
            [state, report] = model.parentModel.updateAfterConvergence(varargin{:});
        end
        
        function dt = getMaximumTimestep(model, state, state0, dt, drivingForces)
            dt = model.parentModel.getMaximumTimestep(state, state0, dt, drivingForces);
        end
        
        function scaling = getScalingFactorsCPR(model, problem, names, solver)
            scaling = model.parentModel.getScalingFactorsCPR(problem, names, solver);
        end
        
        function checkStateFunctionDependencies(model)
            model.parentModel.checkStateFunctionDependencies();
        end
        
        function forces = validateDrivingForces(model, forces, varargin)
            forces = model.parentModel.validateDrivingForces(forces, varargin{:});
        end
        
        function groupings = getStateFunctionGroupings(model)
            groupings = model.parentModel.getStateFunctionGroupings();
        end
        
        function [values, tolerances, names] = getConvergenceValues(model, varargin)
            [values, tolerances, names] = model.parentModel.getConvergenceValues(varargin{:});
        end

        function rmodel = getReservoirModel(model)
            rmodel = model.parentModel;
            while ~isa(rmodel, 'ReservoirModel')
                rmodel = rmodel.parentModel;
            end
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
