classdef NonLinearSolver < handle
%Generalized Newton-like nonlinear solver
%
% SYNOPSIS:
%   solver = NonLinearSolver()
%
%   solver = NonLinearSolver('maxIterations', 5)
%
% DESCRIPTION:
%   The NonLinearSolver class is a general non-linear solver based on
%   Newton's method. It is capable of timestep selection and cutting based
%   on convergence rates and can be extended via subclassing or modular
%   linear solvers and timestep classes.
%
%   Convergence is handled by the PhysicalModel class. The NonLinearSolver
%   simply responds based on what the model reports in terms of convergence
%   to ensure some level of encapsulation.
%
% REQUIRED PARAMETERS:
%   None.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   Documented in methods section.
%
% RETURNS:
%   A NonLinearSolver class instance ready for use.
%
% SEE ALSO:
%   simulateScheduleAD, LinearSolverAD, SimpleTimeStepSelector

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
        % The max number of iterations during a ministep.
        maxIterations
        % The minimum number of solves during a ministep.
        minIterations
        % The maximum number of times the timestep can be halved before it
        % is counted as a failed run
        maxTimestepCuts
        % The solver used to solve the linearized problems during the
        % simulation.
        LinearSolver
        % Verbose flag used to get extra output during simulation.
        verbose
        % Identifier for the nonlinear solver
        identifier
        % Subclass of SimpleTimeStepSelector used to select timesteps
        % during simulation. By default no dynamic timestepping will be
        % enabled.
        timeStepSelector

        % Boolean indicating if Newton increments should be relaxed.
        useRelaxation
        % Relaxation parameter between 0 and 1. This is modified
        % dynamically if useRelaxation is on, and should in general not be
        % modified unless you know what you are doing.
        relaxationParameter
        % Either 'dampen', 'sor' or 'none'
        % For dampen, where w = relaxationParameter.
        %       x_new = x_old + dx*w
        % For successive over-relaxation (SOR)
        %       x_new = x_old + dx*w + dx_prev*(1-w)
        relaxationType
        % Relaxation is reduced by this when stagnation occurs
        relaxationIncrement
        % Lowest possible relaxation factor
        minRelaxation
        % Largest possible relaxation factor
        maxRelaxation
        
        useLinesearch
        alwaysUseLinesearch
        linesearchReductionFactor
        linesearchDecreaseFactor
        linesearchMaxIterations
        linesearchConvergenceNames
        linesearchResidualScaling
        linesearchReductionFn
        
        convergenceIssues

        % Abort a timestep if no reduction is residual is happening.
        enforceResidualDecrease
        % Stagnation tolerance - used in relaxation to determine of a
        % residual value is no longer decreasing
        stagnateTol

        % If error on failure is not enabled, the solver will return even
        % though it did not converge. May be useful for debugging. Results
        % should not be relied upon if this is enabled.
        errorOnFailure
        % If errorOnFailure is disabled, the solver will continue after a
        % failed timestep, treating it as a simply non-converged result
        % with the maximum number of iterations
        continueOnFailure
    end
    
    properties (Access=private)
        % Internal bookkeeping.
        previousIncrement
        previousStepReport
    end

    methods
        function solver = NonLinearSolver(varargin)
            solver.maxIterations   = 25;
            solver.minIterations   = 1;
            solver.identifier      = '';
            solver.verbose         = mrstVerbose();
            solver.maxTimestepCuts = 6;
            solver.LinearSolver    = [];

            solver.relaxationParameter = 1;
            solver.relaxationType = 'dampen';
            solver.useRelaxation = false;
            solver.relaxationIncrement = 0.1;
            solver.minRelaxation = 0.5;
            solver.maxRelaxation = 1;
            
            solver.useLinesearch = false;
            solver.alwaysUseLinesearch = false;
            solver.linesearchReductionFactor = 1/2;
            solver.linesearchDecreaseFactor = 1;
            solver.linesearchMaxIterations = 10;
            solver.linesearchConvergenceNames = {};
            solver.linesearchResidualScaling = 1;
            solver.linesearchReductionFn = [];
            
            solver.enforceResidualDecrease = false;
            solver.stagnateTol = 1e-2;
            solver.convergenceIssues = false;

            solver.errorOnFailure = true;
            solver.continueOnFailure = false;

            solver = merge_options(solver, varargin{:});

            if isempty(solver.LinearSolver)
                solver.LinearSolver = BackslashSolverAD();
            end

            if isempty(solver.timeStepSelector)
                solver.timeStepSelector = SimpleTimeStepSelector();
            end
        end

        function [state, report, ministates] = solveTimestep(solver, state0, dT, model, varargin)
            % Solve a timestep for a non-linear system using one or more substeps
            % REQUIRED PARAMETERS:
            %   state0    - State at the beginning of the timestep
            %   dT        - Timestep size. The solver will move forwards
            %               either as a single step or multiple substeps
            %               depending on convergence rates and sub timestep
            %               selection.
            %   model     - Model inheriting from PhysicalModel with a
            %               valid implementation of the "stepFunction"
            %               member function.
            %
            % OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
            %   'W'       - Wells for the timestep. (struct)
            %   'bc'      - Boundary conditions for the problem (struct).
            %   'src'     - Source terms for the timestep (struct).
            %   
            %   NOTE: Wells, boundary conditions and source terms are the
            %         standard types of external forces in MRST. However,
            %         the model input determines which of these are
            %         actually implemented for that specific step function.
            %         Not all combinations are meaningful for all models.
            %
            %         Some models may implement other types of external
            %         forces that have other names, specified in the
            %         model's "getValidDrivingForces" method.
            %
            % RETURNS:
            %  state      - Problem state after timestep, i.e. if state0
            %               held pressure, saturations, ... at T_0, state
            %               now holds the same values at T_0 + dT.
            %  report     - Report struct, containing some standard
            %               information (iteration count, convergence
            %               status etc) in addition to any reports the
            %               stepFunction contains.
            %  ministates - Cell array containing all ministeps used to get
            %               to T = T_0 + dt. If the solver decided to take
            %               a single step and was successful, this will
            %               just be {state}.
            % SEE ALSO:
            %   PhysicalModel

            opt = struct('initialGuess', state0);

            % Get default driving forces for model
            drivingForces = model.getValidDrivingForces();
            % Add optional control ID for checking when forces are changing
            drivingForces.controlId = nan;

            [opt, forcesArg] = merge_options(opt, varargin{:});
            % Merge in forces as varargin
            drivingForces = merge_options(drivingForces, forcesArg{:});
            model = model.validateModel(drivingForces);

            assert(dT >= 0, [solver.getId(), 'Negative timestep detected.']);

            converged = false;
            done = false;

            dt = dT;

            % Number of nonlinear iterations total
            itCount = 0;
            % Number of ministeps due to cutting
            cuttingCount = 0;
            % Number of steps
            stepCount = 0;
            % Number of accepted steps
            acceptCount = 0;

            t_local = 0;

            isFinalMinistep = false;
            state0_inner = state0;
            % Previous state for a given timestep
            state_prev = [];
            % Timestep that led from state_prev to current state0_inner
            dt_prev = nan;

            wantMinistates = nargout > 2;
            [reports, ministates] = deal(cell(min(2^solver.maxTimestepCuts, 128), 1));

            state = opt.initialGuess;

            % Let the step selector know that we are at start of timestep
            % and what the current driving forces are
            stepsel = solver.timeStepSelector;
            stepsel.newControlStep(drivingForces);

            dtMin = dT/(2^solver.maxTimestepCuts);
            while ~done
                dt_sel = stepsel.pickTimestep(dt_prev, dt, model, solver, state_prev, state0_inner);
                if t_local + dt_sel >= dT
                    % Ensure that we hit report time
                    isFinalMinistep = true;
                    dt = dT - t_local;
                else
                    dt = dt_sel;
                end
                if solver.verbose && dt < dT
                    fprintf('%sSolving ministep : %s (%1.2f %% of control step, control step currently %1.2f %% complete)\n',...
                            solver.getId(), formatTimeRange(dt), dt / dT * 100, t_local / dT * 100)
                end
                [state, failure, tmp] = ...
                    solveMinistep(solver, model, state, state0_inner, dt, drivingForces);

                % Store timestep info
                converged = tmp.Converged;
                its = tmp.Iterations;
                tmp.LocalTime = t_local + dt;
                reports{end+1} = tmp; %#ok
                clear tmp
                if isFinalMinistep && dt/dt_sel > 0.9
                    % Avoid storing ministeps that are just due to cutting
                    % at the end of the control step
                    stepsel.storeTimestep(reports{end});
                end

                % Keep total itcount so we know how much time we are
                % wasting
                itCount = itCount + its;
                stepCount = stepCount + 1;
                if converged
                    t_local = t_local + dt;
                    dt_prev = dt;
                    state_prev = state0_inner;
                    state0_inner = state;
                    acceptCount = acceptCount + 1;

                    if wantMinistates
                        % Output each substep
                        nm = numel(ministates);
                        if nm < acceptCount
                            tmp = cell(nm*2, 1);
                            tmp(1:nm) = ministates;
                            ministates = tmp;
                            clear tmp
                        end
                        ministates{acceptCount} = state;
                    end
                else
                    % Model did not converge, we are in some kind of
                    % trouble.
                    stopNow = dt <= dtMin || failure;
                    
                    if ~(stopNow && solver.continueOnFailure)
                        if acceptCount == 0
                            % We are still at the beginning, and we must honor
                            % the initial guess for current state
                            state = opt.initialGuess;
                        else
                            % Otherwise we reset the initial guess to the
                            % previous step, as we are between ministeps.
                            state = state0_inner;
                        end
                    end
                    msg = [solver.getId(), 'Did not find a solution: '];
                    if failure
                        % Failure means something is seriously wrong,
                        % and we should abort the entire control step
                        % immediately. The last report should include a
                        % FailureMsg field that tells the user what
                        % went wrong.
                        msg = [msg, 'Model step resulted in failure state. Reason: ', ...
                               reports{end}.NonlinearReport{end}.FailureMsg];
                    else
                        msg = [msg, 'Maximum number of substeps stopped timestep reduction'];
                    end

                    if stopNow
                        if solver.errorOnFailure
                            error(msg);
                        else
                            if solver.verbose >= 0
                                warning(msg);
                            end
                            converged = false;
                            if ~(solver.continueOnFailure && failure)
                                % Unless this was a failure and a special
                                % option was set, we are all done here.
                                break;
                            end
                        end
                    else
                        if solver.verbose >= 0
                            % Beat timestep with a hammer
                            warning([solver.getId(), 'Solver did not converge, cutting timestep'])
                        end
                        cuttingCount = cuttingCount + 1;
                        dt = dt/2;
                    end
                    isFinalMinistep = false;
                end
                done = isFinalMinistep && converged;
            end

            if acceptCount ~= 1, pl_mini = 's'; else pl_mini = ''; end
            if itCount     ~= 1, pl_it   = 's'; else pl_it   = ''; end

            dispif(solver.verbose > 0, ...
                   [solver.getId(), ...
                    'Solved timestep with %d accepted ministep%s', ...
                    ' (%d rejected, %d total iteration%s)\n'], ...
                   acceptCount, pl_mini, stepCount - acceptCount, ...
                   itCount, pl_it);

            % Truncate reports from step functions
            reports = reports(~cellfun(@isempty, reports));
            report = struct('Iterations',           itCount,...
                            'Converged',            converged,...
                            'MinistepCuttingCount', cuttingCount);
            % Add seperately because struct constructor interprets cell
            % arrays as repeated structs.
            report.StepReports = reports;
            if wantMinistates
                ministates = ministates(~cellfun(@isempty, ministates));
            end
        end

        function [state, failure, report] = solveMinistep(solver, model, state, state0, dt, drivingForces)
            % Attempt to solve a single mini timestep while trying to avoid
            % stagnation or oscillating residuals.
            omega0 = solver.relaxationParameter;
            solver.convergenceIssues = false;
            if model.stepFunctionIsLinear
                maxIts = 0;
            else
                maxIts = solver.maxIterations;
            end
            reports = cell(maxIts, 1);
            for i = 1:(maxIts + 1)
                % If we are past maximum number of iterations, step function will
                % just check convergence and return
                [state, stepReport] = ...
                    model.stepFunction(state, state0, dt, drivingForces, ...
                    solver.LinearSolver, solver, ...
                    i);
                converged  = stepReport.Converged;
                failure    = stepReport.Failure;
                reports{i} = stepReport;
                if converged
                    break
                end
                if failure
                    break
                end

                if i > 1 && solver.enforceResidualDecrease
                    if all(stepReport.Residuals >= prev_best)
                        % We are not seeing reduction, but rather increase in the
                        % residuals. Break and let the solver decide to either
                        % abort or cut the timestep.
                        break;
                    end
                end
                prev_best = stepReport.Residuals;
                
                if solver.useRelaxation || solver.useLinesearch
                    % Store residual history during nonlinear loop to detect
                    % stagnation or oscillations in residuals.
                    if i == 1
                        res = nan(maxIts + 1, numel(stepReport.Residuals));
                    end
                    res(i, :) = stepReport.Residuals;

                    isOk = stepReport.ResidualsConverged;
                    isOscillating = solver.checkForOscillations(res, i);
                    isStagnated = solver.checkForStagnation(res, i);
                    % We will use relaxations if all non-converged residuals are
                    % either stagnating or oscillating.
                    bad = (isOscillating | isStagnated) | isOk;
                    relax = all(bad) && ~all(isOk);
                    if relax
                        if solver.verbose > 0 && ~solver.convergenceIssues
                            fprintf('Convergence issues detected.');
                            if solver.useLinesearch
                                fprintf(' Activating line search.\n');
                            else
                                fprintf(' Activating relaxation.\n');
                            end
                        end
                        solver.convergenceIssues = true;
                        solver.relaxationParameter = max(solver.relaxationParameter - solver.relaxationIncrement, solver.minRelaxation);
                    else
                        solver.relaxationParameter = min(solver.relaxationParameter + solver.relaxationIncrement, solver.maxRelaxation);
                    end
                end
            end
            % If we converged, the last step did not solve anything
            its = i - converged;
            reports = reports(~cellfun(@isempty, reports));
            solver.relaxationParameter = omega0;
            solver.convergenceIssues = false;
            if converged
                [state, r] = model.updateAfterConvergence(state0, state, dt, drivingForces);
                reports{end}.FinalUpdate = r;
            end
            report = struct('NonlinearReport', {reports}, ...
                            'Converged',       converged, ...
                            'Timestep',        dt, ...
                            'Iterations',      its);
        end

        function [dx, report] = stabilizeNewtonIncrements(solver, model, problem, dx)
            % Attempt to stabilize newton increment by changing the values
            % of the increments.
            dx_prev = solver.previousIncrement;
            w = solver.relaxationParameter;
            report = struct('relaxationParameter', w);
            if w < 1
                switch(lower(solver.relaxationType))
                  case 'dampen'
                    for i = 1:numel(dx)
                        dx{i} = dx{i}*w;
                    end
                  case 'sor'
                    if isempty(dx_prev)
                        return
                    end
                    for i = 1:numel(dx)
                        dx{i} = dx{i}*w + (1-w)*dx_prev{i};
                    end
                  case 'none'

                  otherwise
                    error('Unknown relaxationType: Valid options are ''dampen'', ''none'' or ''sor''');
                end
            end
            solver.previousIncrement = dx;
            
        end

        function [stateNext, updateReport, lineReport] = applyLinesearch(solver, model, state0, state, problem0, dx, drivingForces, varargin)
            assert(solver.linesearchReductionFactor < 1 & solver.linesearchReductionFactor > 0, ...
                    'NonLinearSolver.linesearchReductionFactor must be less than unity and positive.');
            iteration = problem0.iterationNo;
            dt = problem0.dt;
            % Function handle for assembling system equations
            assemble = @(state)  model.getEquations(state0, state, dt, drivingForces, ...
                       'ResOnly', true, ...
                       'iteration', iteration, ...
                       varargin{:});
            % Function for computing updated values with a given delta
            update = @(dx) model.updateState(state, problem0, dx, drivingForces);
            factor = solver.linesearchDecreaseFactor;
            converged = false;
            
            % Check convergence of previous iteration. This is the value to
            % beat, i.e. we want a reduction in the residual with respect
            % to this value.
            [ok, val0, names] = model.checkConvergence(problem0);
            activeNames = getActiveNames(solver, names);
            vBest = linesearchApplyUpdate(solver, val0, ok, activeNames);
            
            for its = 1:solver.linesearchMaxIterations
                [stateNext, updateReport] = update(dx);
                problem = assemble(stateNext);

                [ok, val] = model.checkConvergence(problem);
                v = linesearchApplyUpdate(solver, val, ok, activeNames);
                if all(ok) || (any(v < vBest*factor) && all(v <= vBest))
                    dispif(solver.verbose, 'Linesearch reduction successful after %d steps!\n', its)
                    converged = true;
                    break
                end
                % Multiply the increments with the reduction factor
                dx = cellfun(@(x) x.*solver.linesearchReductionFactor, dx, 'UniformOutput', false);
            end
            lineReport = struct('Iterations', its, ...
                                'Converged', converged);
            dispif(solver.verbose && ~converged, 'Linesearch was unable to reduce residual.\n');
        end

        function isOscillating = checkForOscillations(solver, res, index) %#ok
        % Check if residuals are oscillating. They are oscillating of
        % the ratio of forward and backwards differences for a specific
        % residual is negative.
            if index < 3
                isOscillating = false(1, size(res, 2));
                return
            end

            old = res(index - 2, :);
            mid = res(index - 1, :);
            next = res(index,    :);
            dfdb = (next - mid)./(mid - old);
            isOscillating = dfdb < 0;
        end

        function isStagnated = checkForStagnation(solver, res, index)
        % Check if residuals have stagnated. Residuals are flagged as
        % stagnating if the relative change is smaller than
        % the tolerance (in absolute value).
            if index < 2
                isStagnated = false(1, size(res, 2));
                return
            end
            prev = res(index - 1, :);
            next = res(index,     :);

            isStagnated = abs(next - prev)./prev < solver.stagnateTol;
        end

        function str = getId(solver)
            if isempty(solver.identifier)
                str = '';
            else
                str = [solver.identifier, ': '];
            end
        end
    end
end

function v = linesearchApplyUpdate(solver, v, ok, active)
    v = v./solver.linesearchResidualScaling;
    v = v(active);
    ok = ok(active);

    if ~isempty(solver.linesearchReductionFn)
        % Apply function
        v = solver.linesearchReductionFn(v);
    else
        % Set converged value to zero
        v = v.*~ok;
    end
end

function activeNames = getActiveNames(solver, names)
    if isempty(solver.linesearchConvergenceNames)
        activeNames = true(size(names));
    else
        activeNames = ismember(names, solver.linesearchConvergenceNames);
    end
end
