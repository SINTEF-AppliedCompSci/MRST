classdef PhysicalModel
%Base class for physical models
%
% SYNOPSIS:
%   model = PhysicalModel(G)
%
% DESCRIPTION:
%   Base class for implementing physical models for use with automatic
%   differentiation. This class cannot be used directly.
%
%   A physical model consists of a set of discrete operators that can be
%   used to define the model equations and a nonlinear tolerance that
%   defines how close the values must be to zero before the equations can
%   be considered to be fulfilled. In most cases, the operators are defined
%   over a grid, which is an optional property in this class. In addition,
%   the class contains a flag informing if the model equations are linear,
%   and a flag determining verbosity of class functions.
%  
%   The class contains member functions for:
%     - evaluating residual equations and Jacobians
%     - querying and setting individual variables in the physical state
%     - executing a single nonlinear step (i.e., a linear solve with a
%       possible subsequent stabilization step), verifying convergence, and
%       reporting the status of the step
%     - verifying the model, associated physical states, or individual
%       physical properties
%   as well as a number of utility functions for updating the physical
%   state with increments from the linear solve, etc. See the
%   implementation of the class for more details.
%
% REQUIRED PARAMETERS:
%
%   G  - Simulation grid.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   See class properties.
%
% RETURNS:
%   Class instance.
%
% SEE ALSO:
%   ReservoirModel, ThreePhaseBlackOilModel, TwoPhaseOilWaterModel

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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

properties
    % Operators used for construction of systems
    operators
    % Inf norm tolerance for nonlinear iterations
    nonlinearTolerance
    % Grid. Can be empty.
    G
    % Verbosity from model routines
    verbose
    % Model step function is guaranteed to converge in a single step.
    % Do not enable this unless you are very certain that it is the
    % case!
    stepFunctionIsLinear
end

methods
    function model = PhysicalModel(G, varargin)
        model.nonlinearTolerance = 1e-6;
        model.verbose = mrstVerbose();
        model = merge_options(model, varargin{:});
        model.G = G;

        model.stepFunctionIsLinear = false;
    end

    % --------------------------------------------------------------------%
    function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin) %#ok
        % Get the set of linearized equations governing the system; this
        % should be an instance of the LinearizedProblem class containing
        % the residual equations + Jacobians etc for the model. We also
        % return the state, because the equation setup can in some special
        % cases modify the state.
        error('Base class not meant for direct use')
    end

    % --------------------------------------------------------------------%
    function [problem, state] = getAdjointEquations(model, state0, state, dt, drivingForces, varargin)
        % Function to get equation when using adjoint to calculate
        % gradients. This make it possible to use different equations to
        % calculate the solution in the forward mode forexample if
        % equations are solved explicitely like for hysteretic models.
        % it is assumed that the solution of the system in forward fot the
        % two diffent equations are equal i.e problem.val == 0 .
        [problem, state] = model.getEquations(state0, state, dt, drivingForces, varargin{:});
    end

    % --------------------------------------------------------------------%
    function state = validateState(model, state) %#ok
        % Validate the state for use with the model. Should check that
        % required fields are present and of the right dimensions. If
        % missing fields can be assigned default values, state is return
        % with the required fields added. If reasonable default values
        % cannot be assigned, a descriptive error should be thrown telling
        % the user what is missing or wrong (and ideally how to fix it).

        % Any state is valid for base class
        return
    end

    % --------------------------------------------------------------------%
    function model = validateModel(model, varargin)
        % Validate that a model is suitable for simulation. If the missing
        % or inconsistent parameters can be fixed automatically, an updated
        % model will be returned. Otherwise, an error should occur.
        %
        % Second input may be the forces struct argument. This function
        % should NOT require forces arg to run, however.

        % Base class is always suitable
        return
    end

    % --------------------------------------------------------------------%
    function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
    % Update state based on Newton increments
        for i = 1:numel(problem.primaryVariables);
             p = problem.primaryVariables{i};
             % Update the state
             state = model.updateStateFromIncrement(state, dx, problem, p);
        end
        report = [];
    end

    % --------------------------------------------------------------------%
    function [model, state] = updateForChangedControls(model, state, drivingForces) %#ok
        % Whenever controls change, this function should ensure that both
        % model and state are up to date with the present set of driving
        % forces.
    end
    % --------------------------------------------------------------------%
    function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces) %#ok
        % Update state based on nonlinear increment after timestep has
        % converged. Defaults to doing nothing since not all models
        % require this.
        report = [];
    end

    % --------------------------------------------------------------------%
    function [convergence, values, names] = checkConvergence(model, problem, n)
        % Check and report convergence based on residual tolerances
        if nargin == 2
            n = inf;
        end

        values = norm(problem, n);
        convergence = values < model.nonlinearTolerance;
        names = strcat(problem.equationNames, ' (', problem.types, ')');
    end

    % --------------------------------------------------------------------%
    function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin)
        % Execute a single linearized timestep. Models with linear flow
        % equations are given special treatment.

        onlyCheckConvergence = iteration > nonlinsolver.maxIterations;

        [problem, state] = model.getEquations(state0, state, dt, drivingForces, ...
                                   'ResOnly', onlyCheckConvergence, ...
                                   'iteration', iteration, ...
                                   varargin{:});
        problem.iterationNo = iteration;
        problem.drivingForces = drivingForces;
        
        [convergence, values, resnames] = model.checkConvergence(problem);
        
        % Minimum number of iterations can be prescribed, i.e., we
        % always want at least one set of updates regardless of
        % convergence criterion.
        doneMinIts = iteration > nonlinsolver.minIterations;

        % Defaults
        failureMsg = '';
        failure = false;
        [linearReport, updateReport, stabilizeReport] = deal(struct());
        if (~(all(convergence) && doneMinIts) && ~onlyCheckConvergence)
            % Get increments for Newton solver
            [dx, ~, linearReport] = linsolver.solveLinearProblem(problem, model);
            if any(cellfun(@(d) ~all(isfinite(d)), dx))
                failure = true;
                failureMsg = 'Linear solver produced non-finite values.';
            end
            % Let the nonlinear solver decide what to do with the
            % increments to get the best convergence
            [dx, stabilizeReport] = nonlinsolver.stabilizeNewtonIncrements(model, problem, dx);

            if (nonlinsolver.useLinesearch && nonlinsolver.convergenceIssues) || ...
                nonlinsolver.alwaysUseLinesearch
                [state, updateReport, stabilizeReport.linesearch] = nonlinsolver.applyLinesearch(model, state0, state, problem, dx, drivingForces, varargin{:});
            else
                % Finally update the state. The physical model knows which
                % properties are actually physically reasonable.
                [state, updateReport] = model.updateState(state, problem, dx, drivingForces);
            end
        end
        isConverged = (all(convergence) && doneMinIts) || model.stepFunctionIsLinear;
        
        % If step function is linear, we need to call a residual-only
        % equation assembly to ensure that indirect/derived quantities are
        % set with the updated values (fluxes, mobilities and so on).
        if model.stepFunctionIsLinear
            [~, state] = model.getEquations(state0, state, dt, drivingForces, ...
                                   'ResOnly', true, ...
                                   'iteration', iteration+1, ...
                                   varargin{:});
        end
        if model.verbose
            printConvergenceReport(resnames, values, convergence, iteration);
        end
        report = model.makeStepReport(...
                        'LinearSolver', linearReport, ...
                        'UpdateState',  updateReport, ...
                        'Failure',      failure, ...
                        'FailureMsg',   failureMsg, ...
                        'Converged',    isConverged, ...
                        'Residuals',    values, ...
                        'StabilizeReport', stabilizeReport,...
                        'ResidualsConverged', convergence);
    end

    % --------------------------------------------------------------------%
    function report = makeStepReport(model, varargin) %#ok
        % Get the standardized step report that all models produce.
        report = struct('LinearSolver', [], ...
                        'UpdateState',  [], ...
                        'Failure',      false, ...
                        'FailureMsg',   '', ...
                        'Converged',    false, ...
                        'FinalUpdate',  [],...
                        'Residuals',    [],...
                        'StabilizeReport', [], ...
                        'ResidualsConverged', []);
        report = merge_options(report, varargin{:});
    end

    % --------------------------------------------------------------------%
    function [gradient, result, report] = solveAdjoint(model, solver, getState,...
                                getObjective, schedule, gradient, itNo)
        % Solve the adjoint equations for a given step.
        validforces = model.getValidDrivingForces();
        dt_steps = schedule.step.val;

        current = getState(itNo);
        before    = getState(itNo - 1);        
        dt = dt_steps(itNo);

        lookupCtrl = @(step) schedule.control(schedule.step.control(step));
        % get forces and merge with valid forces
        forces = model.getDrivingForces(lookupCtrl(itNo));
        forces = merge_options(validforces, forces{:});
        model = model.validateModel(forces);
        
        % Initial state typically lacks wellSol-field, so add if needed
        if itNo == 1
            before = model.validateState(before);
        end
        
        % We get the forward equations via the reverseMode flag. This is
        % slightly hacky (we should use the forward equations instead), but
        % we assume that the reverseMode flag takes care of this for us.
        % This slightly messy setup is made to support hysteresis models,
        % where the forward and backwards equations are very different.
        problem = model.getAdjointEquations(before, current, dt, forces, ...
                                    'reverseMode', false, 'iteration', inf);

        if itNo < numel(dt_steps)
            after    = getState(itNo + 1);
            dt_next = dt_steps(itNo + 1);
            % get forces and merge with valid forces
            forces_p = model.getDrivingForces(lookupCtrl(itNo + 1));
            forces_p = merge_options(validforces, forces_p{:});
            problem_p = model.getAdjointEquations(current, after, dt_next, forces_p,...
                                'iteration', inf, 'reverseMode', true);
        else
            problem_p = [];
        end
        [gradient, result, rep] = solver.solveAdjointProblem(problem_p,...
                                    problem, gradient, getObjective(itNo), model);
        report = struct();
        report.Types = problem.types;
        report.LinearSolverReport = rep;
    end

    % --------------------------------------------------------------------%
    function [fn, index] = getVariableField(model, name)                   %#ok
        % Get the index/name mapping for the model (such as where
        % pressure or water saturation is located in state). This
        % always result in an error, as this model knows of no variables.
        %
        % Given a name, this function produces the fieldname and the
        % index in the struct that can be used to get the same info.
        [fn, index] = deal([]);

        if isempty(index)
            error('PhysicalModel:UnknownVariable', ...
                ['State variable ''', name, ''' is not known to this model']);
        end
    end

    % --------------------------------------------------------------------%
    function p = getProp(model, state, name)
        % Get a property based on the name. Uses getVariableField to
        % determine how to obtain the data for the name, e.g.,
        %    p = model.getProp(state, 'pressure');
        
        [fn, index] = model.getVariableField(name);
        p = state.(fn)(:, index);
    end

    % --------------------------------------------------------------------%
    function varargout = getProps(model, state, varargin)
        % Get multiple properties based on their name(s). Multiple
        % names can be sent in as variable arguments, e.g.,
        %    [p, s] = model.getProps(state, 'pressure', 's');

        varargout = cellfun(@(x) model.getProp(state, x), ...
                            varargin, 'UniformOutput', false);
    end

    % --------------------------------------------------------------------%
    function state = incrementProp(model, state, name, increment)
        % Increment property based on the name for the field. The
        % returned state contains incremented values. E.g., 
        %   state = struct('pressure', 0);
        %   state = model.incrementProp(state, 'pressure', 1);
        % will set state.pressure to 1 if the variable pressure is known to
        % the model. Otherwise we will get an error.
        
        [fn, index] = model.getVariableField(name);
        p = state.(fn)(:, index)  + increment;
        state.(fn)(:, index) = p;
    end

    % --------------------------------------------------------------------%
    function state = setProp(model, state, name, value)
        % Set property to given value based on name. E.g.,
        %   state = struct('pressure', 0);
        %   state = model.setProp(state, 'pressure', 5);
        % will set state.pressure to 5, unless it is not a valid field for
        % the specific model.
        
        [fn, index] = model.getVariableField(name);
        state.(fn)(:, index) = value;
    end

    % --------------------------------------------------------------------%
    function dv = getIncrement(model, dx, problem, name)                   %#ok
        % Find increment in linearized problem with given name, or
        % output zero if not found. A linearized problem can give
        % updates to multiple variables and this makes it easier to get
        % those values without having to know the order they were input
        % into the constructor.

        isVar = problem.indexOfPrimaryVariable(name);
        if any(isVar)
            dv = dx{isVar};
        else
            dv = 0;
        end
    end

    % --------------------------------------------------------------------%
    function [state, val, val0] = updateStateFromIncrement(model, state, dx, problem, name, relchangemax, abschangemax)
        % Update a state, with optionally maximum changes (relative and
        % absolute). E.g., 
        %
        %   state = struct('pressure', 10);
        %   state = model.updateStateFromIncrement(state, 100, problem, 'pressure')
        %
        % will result in pressure being set to 110. Consider alternatively:
        %
        %   state = model.updateStateFromIncrement(state, 100, problem, 'pressure', .1) 
        %
        % which will set pressure to 11, as setting it directly to 110
        % will violate the maximum relative change of 10% 
        %
        % Relative limits such as these are important when working with
        % tabulated and nonsmooth properties in a Newton-type loop, as the
        % initial updates may be far outside the reasonable region of
        % linearization for a complex problem. On the other hand, limiting
        % the relative updates can delay convergence for smooth problems
        % with analytic properties and will, in particular, prevent zero
        % states from being updated, so use with care.
        
        if iscell(dx)
            % We have cell array of increments, use the problem to
            % determine where we can actually find it.
            dv = model.getIncrement(dx, problem, name);
        else
            % Numerical value, increment directly and do not safety
            % check that this is a part of the model
            dv = dx;
        end

        val0 = model.getProp(state, name);

        [changeRel, changeAbs] = deal(1);
        if nargin > 5
            [~, changeRel] = model.limitUpdateRelative(dv, val0, relchangemax);
        end
        if nargin > 6
            [~, changeAbs] = model.limitUpdateAbsolute(dv, abschangemax);
        end            
        % Limit update by lowest of the relative and absolute limits 
        change = min(changeAbs, changeRel);

        val   = val0 + dv.*repmat(change, 1, size(dv, 2));
        state = model.setProp(state, name, val);
    end
    % --------------------------------------------------------------------%
    function state = capProperty(model, state, name, minvalue, maxvalue)
        % Cap values to min/max values
        v = model.getProp(state, name);
        v = max(minvalue, v);
        if nargin > 4
            v = min(v, maxvalue);
        end
        state = model.setProp(state, name, v);
    end
    
    % --------------------------------------------------------------------%
    function [vararg, control] = getDrivingForces(model, control) %#ok
        % Setup and pass on driving forces. Outputs both a cell array
        % suitable for parsing by merge_options and the control struct.
        vrg = [fieldnames(control), struct2cell(control)];
        vararg = reshape(vrg', 1, []);
    end
    
    % --------------------------------------------------------------------%
    function forces = getValidDrivingForces(model) %#ok
        % Get struct with valid forces for model, with reasonable default
        % values.
        forces = struct();
    end
    
    % --------------------------------------------------------------------%
    function checkProperty(model, state, property, n_el, dim)
        % model    - model to be checked
        % state    - state to be checked
        % property - valid property (according to getVariableField)
        % n_el     - array of same size as dim, indicating how many entries
        %            there should be along each dimension.
        % dim      - dimensions to check.
        if numel(dim) > 1
            assert(numel(n_el) == numel(dim));
            % Recursively check all dimensions 
            for i = 1:numel(dim)
                model.checkProperty(state, property, n_el(i), dim(i));
            end
            return
        end
        fn = model.getVariableField(property);
        assert(isfield(state, fn), ['Field ".', fn, '" missing! ', ...
            property, ' must be supplied for model "', class(model), '"']);

        if dim == 1
            sn = 'rows';
        elseif dim == 2
            sn = 'columns';
        else
            sn = ['dimension ', num2str(dim)];
        end
        n_actual = size(state.(fn), dim);
        assert(n_actual == n_el, ...
            ['Dimension mismatch for ', sn, ' of property "', property, ...
            '" (state.', fn, '): Expected ', sn, ' to have ', num2str(n_el), ...
            ' entries but state had ', num2str(n_actual), ' instead.'])
    end
end

methods (Static)
    % --------------------------------------------------------------------%
    function [dv, change] = limitUpdateRelative(dv, val, maxRelCh)
        % Limit a update by relative limit
        biggestChange = max(abs(dv./val), [], 2);
        change = min(maxRelCh./biggestChange, 1);
        dv = dv.*repmat(change, 1, size(dv, 2));
    end
    % --------------------------------------------------------------------%
    function [dv, change] = limitUpdateAbsolute(dv, maxAbsCh)
        % Limit a update by absolute limit
        biggestChange = max(abs(dv), [], 2);
        change = min(maxAbsCh./biggestChange, 1);
        dv = dv.*repmat(change, 1, size(dv, 2));
    end
    % --------------------------------------------------------------------%
    function [vars, isRemoved] = stripVars(vars, names)
        isRemoved = cellfun(@(x) any(strcmpi(names, x)), vars);
        vars(isRemoved) = [];
    end
end
end

