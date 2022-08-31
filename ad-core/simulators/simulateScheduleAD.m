function [wellSols, states, schedulereport] = ...
      simulateScheduleAD(initState, model, schedule, varargin)
% Run a schedule for a non-linear physical model using an automatic differention
%
% SYNOPSIS:
%   wellSols = simulateScheduleAD(initState, model, schedule)
%
%   [wellSols, state, report]  = simulateScheduleAD(initState, model, schedule)
%
% DESCRIPTION:
%   This function takes in a valid schedule file (see required parameters)
%   and runs a simulation through all timesteps for a given model and
%   initial state.
%
%   simulateScheduleAD is the outer bookkeeping routine used for running
%   simulations with non-trivial control changes. It relies on the model
%   and (non)linear solver classes to do the heavy lifting.
%
% REQUIRED PARAMETERS:
%
%   initState    - Initial reservoir/model state. It should have whatever
%                  fields are associated with the physical model, with
%                  reasonable values. It is the responsibility of the user
%                  to ensure that the state is properly initialized.
%
%   model        - The physical model that determines jacobians/convergence
%                  for the problem. This must be a subclass of the
%                  `PhysicalModel` base class.
%
%   schedule     - Schedule containing fields step and control, defined as
%                  follows:
%                         - `schedule.control` is a struct array containing
%                           fields that the model knows how to process.
%                           Typically, this will be the fields such as `.W` 
%                           for wells or `.bc` for boundary conditions.
%
%                         - `schedule.step` contains two arrays of equal 
%                           size named `val` and `control`. Control is a
%                           index into the `schedule.control` array,
%                           indicating which control is to be used for the
%                           timestep.`schedule.step.val` is the timestep
%                           used for that control step.
%
% OPTIONAL PARAMETERS:
%
%   'Verbose'         - Indicate if extra output is to be printed such as
%                       detailed convergence reports and so on.
%
%   'OutputMinisteps' - The solver may not use timesteps equal to the
%                       control steps depending on problem stiffness and
%                       timestep selection. Enabling this option will make
%                       the solver output the states and reports for all
%                       steps actually taken and not just at the control
%                       steps. See also `convertReportToSchedule` which can
%                       be used to construct a new schedule from these
%                       timesteps.
%
%   'initialGuess'     - A cell array with one entry per control-step. If
%                        provided, the state from this cell array is passed as
%                        initial guess for the `NonLinearSolver`. See:
%                        `NonLinearSolver.solveTimestep`, optional arguments.
%
%   'NonLinearSolver'  - An instance of the `NonLinearSolver` class. Consider
%                        using this if you for example want a special timestep
%                        selection algorithm. See the `NonLinearSolver` class
%                        docs for more information.
%
%   'OutputHandler'    - Output handler class, for example for writing states
%                        to disk during the simulation or in-situ
%                        visualization. See the ResultHandler base class.
%
%   'WellOutputHandler'- Same as 'OutputHandler', but for the well
%                        solutions for the individual report steps. Well
%                        solutions are also stored using OutputHandler, but
%                        using WellOutputHandler is convenient for quickly
%                        loading well solutions only.
%
%   'ReportHandler'    - Same as 'OutputHandler', but for the reports for the
%                        individual report steps.
%
%   'LinearSolver'     - Class instance subclassed from `LinearSolverAD`. Used
%                        to solve linearized problems in the `NonLinearSolver`
%                        class. Note that if you are passing a
%                        `NonLinearSolver`, you might as well put it there.
%
%   'afterStepFn'      - Function handle to an optional function that will be
%                        called after each successful report step in the
%                        schedule. The function should take in the following
%                        input arguments:
%                        - model: The model used in the schedule
%                        - states: A cell array of all states that are
%                          computed, as well as possible empty entries
%                          where the states have not been computed yet.
%                        - reports: A cell array of reports for each step,
%                          with empty entries for steps that have not been
%                          reached yet.
%                        - solver: The NonLinearSolver instance.
%                        - schedule: The current schedule.
%                        - simtime: Array with the time in seconds taken by
%                          the `NonLinearSolver` to compute each step.
%                          Entries not computed will contain zeros.
%
%                        See `getPlotAfterStep` for more information and
%                        `blackoilTutorialPlotHook` for a worked example.
%
%   'processOutputFn'  - Function handle to an optional function that
%                        processes the simulation output (wellSols, states
%                        and reports) before they are stored to state using
%                        the output handlers. Allows for storing only data
%                        of interest, which is useful when dealing with
%                        large models. Changes made to the output by this
%                        function are only applied to the data that is
%                        stored, and will not affect what is passed on to
%                        the next timestep.
%
%   
%   'controlLogicFn'   - Function handle to optional function that will be
%                        called after each step enabling schedule updates to 
%                        be triggered on specified events. Input arguemnts:
%                        - state: The current state
%                        - schedule: The current schedule
%                        - report: Current report
%                        - i: The current report step such that current
%                          control step equals schedule.step.control(i)
%                        The function must have three outputs:
%                        - schedule: Possibly updated schedule
%                        - report: Possibly updated report
%                        - isAltered: Flag indicated whether the schedule
%                          was updated or not
%
%   'checkOperators'   - Flag indicating if we should check that model.G
%                        and model.rock are identical to those used to
%                        define model.operators before the simulation
%                        starts. Defaults to true.
%                       
% RETURNS:
%   wellSols         - Well solution at each control step (or timestep if
%                      'OutputMinisteps' is enabled.)
%
%   states           - State at each control step (or timestep if
%                      'OutputMinisteps' is enabled.)
%
%   schedulereport   - Report for the simulation. Contains detailed info for
%                      the whole schedule run, as well as arrays containing
%                      reports for the whole stack of routines called during
%                      the simulation.
%
% NOTE:
%   For examples of valid models, see `ThreePhaseBlackOilModel` or
%   `TwoPhaseOilWaterModel`.
%
% SEE ALSO:
%   `computeGradientAdjointAD`, `PhysicalModel`

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    assert (isa(model, 'PhysicalModel'), ...
            'The model must be derived from PhysicalModel');
    opt = struct('Verbose',           mrstVerbose(),...
                 'OutputMinisteps',   false, ...
                 'initialGuess',      {{}}, ...
                 'NonLinearSolver',   [], ...
                 'OutputHandler',     [], ...
                 'WellOutputHandler', [], ...
                 'ReportHandler',     [], ...
                 'afterStepFn',       [], ...
                 'controlLogicFn',    [], ...
                 'processOutputFn',   [], ...
                 'restartStep',       1,  ...
                 'outputOffset',      [], ...
                 'LinearSolver',      [], ...
                 'checkOperators',    []);

    opt = merge_options(opt, varargin{:});
    if isempty(opt.checkOperators)
        opt.checkOperators = mrstSettings('get', 'useHash');
    end
    
    %----------------------------------------------------------------------
    tm = [0 ; reshape(cumsum(schedule.step.val), [], 1)];
    restart = opt.restartStep;
    if restart ~= 1
        T0 = sum(schedule.step.val(1:restart-1));
        if isfield(initState, 'time') && abs(initState.time - T0) > 10*eps(sum(schedule.step.val))
            warning('Time mismatch in initial state for restart. Expected %s, got %s.\n', ...
                formatTimeRange(T0), formatTimeRange(initState.time));
        end
        nStep = numel(schedule.step.val);
        assert(numel(restart) == 1 && ...
               restart <= nStep &&...
               restart > 1, ...
        ['Restart step must be an index between 1 and ', num2str(nStep), '.']);
        schedule.step.control = schedule.step.control(restart:end);
        schedule.step.val = schedule.step.val(restart:end);
    end
    if opt.Verbose
       simulation_header(numel(tm)-1, tm(end));
    end

    step_header = create_step_header(opt.Verbose, tm, opt.restartStep);

    solver = opt.NonLinearSolver;
    if isempty(solver)
        solver = NonLinearSolver('linearSolver', opt.LinearSolver);
    else
        assert(isa(solver, 'NonLinearSolver'), ...
               ['Passed nonlinear solver is not an instance ', ...
                'of NonLinearSolver class!'])
        if ~isempty(opt.LinearSolver)
            % We got a nonlinear solver, but we still want to override the
            % actual linear solver passed to the higher level schedule function
            % we're currently in.
            assert (isa(opt.LinearSolver, 'LinearSolverAD'), ...
                   ['Passed linear solver is not an instance ', ...
                    'of LinearSolverAD class!'])

            solver.LinearSolver = opt.LinearSolver;
        end
    end
    % Reset timestep selector in case it was used previously.
    solver.timeStepSelector.reset();

    nSteps = numel(schedule.step.val);
    [wellSols, states, reports] = deal(cell(nSteps, 1));
    wantStates = nargout > 1 || ~isempty(opt.afterStepFn);
    wantReport = nargout > 2 || ~isempty(opt.afterStepFn);

    % Check if model is self-consistent and set up for current BC type
    dispif(opt.Verbose, 'Validating model...\n')
    ctrl = schedule.control(schedule.step.control(1));
    [forces, fstruct] = model.getDrivingForces(ctrl);
    model = model.validateModel(fstruct, opt.checkOperators);
    dispif(opt.Verbose, 'Model has been validated.\n')
    % Check dependencies
    dispif(opt.Verbose, 'Checking state functions and dependencies...\n')
    model.checkStateFunctionDependencies();
    dispif(opt.Verbose, 'All checks ok. Model ready for simulation.\n')
    % Validate schedule
    dispif(opt.Verbose, 'Preparing schedule for simulation...\n')
    schedule = model.validateSchedule(schedule);
    dispif(opt.Verbose, 'All steps ok. Schedule ready for simulation.\n')

    % Check if initial state is reasonable
    dispif(opt.Verbose, 'Validating initial state...\n')
    state = model.validateState(initState);
    dispif(opt.Verbose, 'Initial state ok. Ready to begin simulation.\n')
    failure = false;
    simtime = zeros(nSteps, 1);
    prevControl = nan;
    firstEmptyIx = 1;
    for i = 1:nSteps
        step_header(i);
        state0 = state;
        
        currControl = schedule.step.control(i);
        if prevControl ~= currControl
            [forces, fstruct] = model.getDrivingForces(schedule.control(currControl));
            [model, state0]= model.updateForChangedControls(state, fstruct);
            prevControl = currControl;
        end
        
        if isempty(opt.initialGuess)
            extraArg = {};
        else
            guess = model.validateState(opt.initialGuess{i});
            extraArg = {'initialGuess', guess};
        end
        dt = schedule.step.val(i);
        timer = tic();
        if opt.OutputMinisteps
            [state, report, substates] = solver.solveTimestep(state0, dt, model, ...
                                            forces{:}, 'controlId', currControl, extraArg{:});
            % We have potentially several ministeps desired as output and
            % we need to offset the pointer to the first array accordingly.
            ind = firstEmptyIx:(firstEmptyIx + numel(substates) - 1);
            firstEmptyIx = firstEmptyIx + numel(substates);
        else
            [state, report] = solver.solveTimestep(state0, dt, model,...
                                            forces{:}, 'controlId', currControl, extraArg{:});
            % Single state requested for dt. We set values accordingly.
            substates = {state};
            ind = i;
        end

        t = toc(timer);
        simtime(i) = t;
        % Abort simulation
        if ~report.Converged
            warning('NonLinear:Failure', ...
                   ['Nonlinear solver aborted, ', ...
                    'returning incomplete results']);
            failure = true;
            break;
        end

        % Apply control logic
        if ~isempty(opt.controlLogicFn)
            [schedule, report, isAltered] = opt.controlLogicFn(state, schedule, report, i);
            if isAltered
                prevControl = nan; 
            end
        end
        
        if opt.Verbose
           disp_step_convergence(report.Iterations, t);
        end
        % Get wellSols, if they exist
        wellSols_step = getWellSols(substates);
        % Process output before proceeding
        if ~isempty(opt.processOutputFn)
            [substates, wellSols_step, report] = opt.processOutputFn(substates, wellSols_step, report);
        end
        % Store output in handlers, if configured
        writeOutput(opt.OutputHandler, opt, ind, substates)
        writeOutput(opt.WellOutputHandler, opt, ind, wellSols_step)
        writeOutput(opt.ReportHandler, opt, i, report, false)
        
        % Write to the cell arrays that will be outputs from the function
        wellSols(ind) = wellSols_step;
        if wantStates
            states(ind) = substates;
        end
        if wantReport
            reports{i} = report;
        end
        
        if ~isempty(opt.afterStepFn)
            [model, states, reports, solver, ok] = opt.afterStepFn(model, states, reports, solver, schedule, simtime);
            if ~ok
                warning('Aborting due to external function');
                break
            end
        end        
        
        if report.EarlyStop
            break;
        end
    end
    
    if wantReport
        reports = reports(~cellfun(@isempty, reports));

        schedulereport = struct();
        schedulereport.ControlstepReports = reports;
        schedulereport.ReservoirTime = cumsum(schedule.step.val);
        schedulereport.Converged  = cellfun(@(x) x.Converged, reports);
        schedulereport.Iterations = cellfun(@(x) x.Iterations, reports);
        schedulereport.SimulationTime = simtime;
        schedulereport.Failure = failure;
    end
    
    if report.EarlyStop
        nSteps = numel(reports);
        earlyStopMsg = ' (termination triggered by stopFunction) ';
    else
        earlyStopMsg = '';
    end
    
    fprintf('*** Simulation complete. Solved %d control steps in %s%s ***\n',...
                                  nSteps, formatTimeRange((sum(simtime))), earlyStopMsg);
    
end

%--------------------------------------------------------------------------
function ws = getWellSols(states)
    ns = numel(states);
    ws = cell(ns, 1);
    for i = 1:ns
        if isfield(states{i}, 'wellSol')
            ws{i} = states{i}.wellSol;
        end
    end
end

%--------------------------------------------------------------------------

function writeOutput(handler, opt, pos, values, useOffset)
    if nargin < 5
        useOffset = true;
    end
    offset = opt.outputOffset;
    if isempty(offset) || ~useOffset
        offset = opt.restartStep;
    end
    if ~isempty(handler)
        handler(pos + offset - 1) = values; %#ok This is a handle class instance.
    end
end

%--------------------------------------------------------------------------

function simulation_header(nstep, tot_time)
   starline = @(n) repmat('*', [1, n]);

   fprintf([starline(65), '\n', starline(10), ...
            sprintf(' Starting simulation: %5d steps, %5.0f days ', ...
                    nstep, convertTo(tot_time, day)), ...
            starline(9), '\n', starline(65), '\n']);
end

%--------------------------------------------------------------------------

function header = create_step_header(verbose, tm, rsrt)
   nSteps  = numel(tm) - 1;
   nDigits = floor(log10(nSteps)) + 1;

   if verbose
      % Verbose mode.  Don't align '->' token in report-step range
      nChar = 0;
   else
      % Non-verbose mode.  Do align '->' token in report-step range
      nChar = numel(formatTimeRange(tm(end), 2));
   end
   offset = rsrt-1;
   header = @(i) ...
      fprintf('Solving timestep %0*d/%0*d: %-*s -> %s\n', ...
              nDigits, i + offset, nDigits, nSteps, nChar, ...
              formatTimeRange(tm(i + 0 + offset), 2), ...
              formatTimeRange(tm(i + 1 + offset), 2));
end

%--------------------------------------------------------------------------

function disp_step_convergence(its, cputime)
    if its ~= 1
        pl_it = 's';
    else
        pl_it = '';
    end

   fprintf(['Completed %d iteration%s in %2.2f seconds ', ...
            '(%2.2fs per iteration)\n\n'], ...
            its, pl_it, cputime, cputime/its);
end
