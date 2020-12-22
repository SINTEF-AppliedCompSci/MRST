function [wellSols, states, schedulereport] = ...
      simulateScheduleADMICP(initState, model, schedule, varargin)
% Function to run a schedule for the model using an automatic differention
% 
% This function is extended from a file in The MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/ad-core/simulators/simulateScheduleAD.m 
%
% We refer to that function for a complete commented version of the file.
% In this file we coment on the added lines.

%{
Partial copyright 2009-2020, SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2020, NORCE Norwegian Research Centre AS, Computational 
Geosciences and Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
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
                 'restartStep',       1, ...
                 'LinearSolver',      []);

    opt = merge_options(opt, varargin{:});

    tm = [0 ; reshape(cumsum(schedule.step.val), [], 1)];
    if opt.restartStep ~= 1
        nStep = numel(schedule.step.val);
        assert(numel(opt.restartStep) == 1 && ...
               opt.restartStep <= nStep &&...
               opt.restartStep > 1, ...
        ['Restart step must be an index between 1 and ', ...
                                                     num2str(nStep), '.']);
        schedule.step.control = schedule.step.control(opt.restartStep:end);
        schedule.step.val = schedule.step.val(opt.restartStep:end);
    end
    if opt.Verbose
       simulation_header(numel(tm)-1, tm(end));
    end

    step_header = create_step_header(opt.Verbose, tm, opt.restartStep);

    solver = opt.NonLinearSolver;
    if isempty(solver)
        solver = NonLinearSolver('linearSolver', opt.LinearSolver);
    elseif ~isempty(opt.LinearSolver)
        assert (isa(opt.LinearSolver, 'LinearSolverAD'), ...
               ['Passed linear solver is not an instance ', ...
                'of LinearSolverAD class!'])

        solver.LinearSolver = opt.LinearSolver;
    end
    solver.timeStepSelector.reset();

    nSteps = numel(schedule.step.val);
    [wellSols, states, reports] = deal(cell(nSteps, 1));
    wantStates = nargout > 1;
    wantReport = nargout > 2 || ~isempty(opt.afterStepFn);
    dispif(opt.Verbose, 'Validating model...\n')
    ctrl = schedule.control(schedule.step.control(1));
    [forces, fstruct] = model.getDrivingForces(ctrl);
    model = model.validateModel(fstruct);
    dispif(opt.Verbose, 'Model has been validated.\n')
    dispif(opt.Verbose, 'Checking state functions and dependencies...\n')
    model.checkStateFunctionDependencies();
    dispif(opt.Verbose, 'All checks ok. Model ready for simulation.\n')
    dispif(opt.Verbose, 'Preparing schedule for simulation...\n')
    schedule = model.validateSchedule(schedule);
    dispif(opt.Verbose, 'All steps ok. Schedule ready for simulation.\n')
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
            [forces, fstruct] = ...
                     model.getDrivingForces(schedule.control(currControl));
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
            [state, report, ministeps] = solver.solveTimestep(state0, ...
               dt, model,forces{:}, 'controlId', currControl, extraArg{:});
        else
            [state, report] = solver.solveTimestep(state0, dt, model,...
                         forces{:}, 'controlId', currControl, extraArg{:});
        end

        t = toc(timer);
        simtime(i) = t;

        if ~report.Converged
            warning('NonLinear:Failure', ...
                   ['Nonlinear solver aborted, ', ...
                    'returning incomplete results']);
            failure = true;
            break;
        end

        if ~isempty(opt.controlLogicFn)
            [schedule, report, isAltered] = opt.controlLogicFn(state, ...
                                                      schedule, report, i);
            if isAltered
                prevControl = nan; 
            end
        end
        
        if opt.Verbose
           disp_step_convergence(report.Iterations, t);
        end
        if opt.OutputMinisteps
            ind = firstEmptyIx:(firstEmptyIx + numel(ministeps) - 1);
            states_step = ministeps;
        else
            ind = i;
            states_step = {state};
        end

        wellSols_step = cellfun(@(x) x.wellSol, states_step, ...
                                'UniformOutput', false);

        wellSols(ind) = wellSols_step;
        
        if ~isempty(opt.processOutputFn)
            [states_step, wellSols_step, report] =  ...
                   opt.processOutputFn(states_step, wellSols_step, report);
        end
        
        if ~isempty(opt.OutputHandler)
            opt.OutputHandler{ind + opt.restartStep - 1} = states_step;
        end
        
        if ~isempty(opt.WellOutputHandler)
          opt.WellOutputHandler{ind + opt.restartStep - 1} = wellSols_step;
        end
        
        if ~isempty(opt.ReportHandler)
            opt.ReportHandler{i + opt.restartStep - 1} = report;
        end
        firstEmptyIx = firstEmptyIx + numel(states_step);

        if wantStates || ~isempty(opt.afterStepFn)
            states(ind) = states_step;
        end

        if wantReport
            reports{i} = report;
        end
        
        % Stop the simulation when the porosity is less than a given ptol    
        if any(model.rock.poro-states{i}.c-states{i}.b<model.fluid.ptol)
            % We take the previous solution to assure poro > 0
            statess(1:i-1,1)=states(1:i-1,1); 
            states=statess;
            fprintf('Clogging has been reached in at least one cell\n');
            break 
        end
        
        if ~isempty(opt.afterStepFn)
            [model, states, reports, solver, ok] = opt.afterStepFn(...
                       model, states,  reports, solver, schedule, simtime);
            if ~ok
                warning('Aborting due to external function');
                break
            end
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
    % Print the number of steps
    if i < nSteps
        fprintf('Simulation complete. Solved %d control steps in %s\n', ...
                                   i - 1, formatTimeRange((sum(simtime))));
    else 
        fprintf('Simulation complete. Solved %d control steps in %s\n', ...
                                  nSteps, formatTimeRange((sum(simtime))));
    end
end

function simulation_header(nstep, tot_time)
   starline = @(n) repmat('*', [1, n]);

   fprintf([starline(65), '\n', starline(10), ...
            sprintf(' Starting simulation: %5d steps, %5.0f days ', ...
                    nstep, convertTo(tot_time, day)), ...
            starline(9), '\n', starline(65), '\n']);
end

function header = create_step_header(verbose, tm, rsrt)
   nSteps  = numel(tm) - 1;
   nDigits = floor(log10(nSteps)) + 1;

   if verbose
      nChar = 0;
   else
      nChar = numel(formatTimeRange(tm(end), 2));
   end
   offset = rsrt-1;
   header = @(i) ...
      fprintf('Solving timestep %0*d/%0*d: %-*s -> %s\n', ...
              nDigits, i + offset, nDigits, nSteps, nChar, ...
              formatTimeRange(tm(i + 0 + offset), 2), ...
              formatTimeRange(tm(i + 1 + offset), 2));
end

function disp_step_convergence(its, cputime)
   if its ~= 1, pl_it = 's'; else pl_it = ''; end

   fprintf(['Completed %d iteration%s in %2.2f seconds ', ...
            '(%2.2fs per iteration)\n\n'], ...
            its, pl_it, cputime, cputime/its);
end
