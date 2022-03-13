function [ok, status] = simulatePackedProblem(problems, varargin)
%Simulate one or more packed simulation problems
%
% SYNOPSIS:
%   [ok, status] = simulatePackedProblem(problems)
%   [ok, status] = simulatePackedProblem(problems, 'continueOnError',false)
%
% REQUIRED PARAMETERS:
%   problems - Either a single packed problem from 'packSimulationProblem'
%              or multiple as a cell array, e.g. {problem1, problem2, ...}
%
% OPTIONAL PARAMETERS:
%   checkTooMany - Check if there exists more output than report steps.
%                  This can indicate re-definition of a simulation case and
%                  may lead to malformed output.
%
%   continueOnError - If true, errors occurring during simulation will not
%                     result in an error. Useful when trying to simulate
%                     multiple problems in sequence.
%
%   restartStep  - Explicitly specify the restart step, potentially
%                  overwriting any previously simulated results from that
%                  step on.
%
%   plot         - Plot everything using plotPackedProblem after simulation
%                  is done.
%
% RETURNS:
%   ok     - Flag for each case indicating a successful run.
%   status - Structs containing some information about each case (possible
%            error message, time spent, etc).
% EXAMPLE:
%   demoPackedProblems
%
% SEE ALSO:
%   packSimulationProblem, getPackedSimulatorOutput

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

    if isstruct(problems)
        problems = {problems};
    end
    
    assert(iscell(problems));
    
    np = numel(problems);
    opt = struct('checkTooMany',    true, ...
                 'continueOnError', np > 1, ...
 		         'plot',            false, ...
                 'restartStep',     nan);
    opt = merge_options(opt, varargin{:});

    ok = true(np, 1);
    status = repmat(struct('Message', '', 'WallTime', [], 'RestartIndex', nan), np, 1);
    for i = 1:np
        problem = problems{i};
        % Unpack base case
        nls = problem.SimulatorSetup.NonLinearSolver;
        doMinisteps = problem.SimulatorSetup.OutputMinisteps;
        
        state_handler = problem.OutputHandlers.states;
        wellSol_handler = problem.OutputHandlers.wellSols;
        report_handler = problem.OutputHandlers.reports;

        ok = true;

        firstLine = sprintf(' Case "%s" (%s)',...
                problem.BaseName, problem.Name);
        secondLine = sprintf(' Description: "%s"',...
                problem.Description);
        
        n1 = numel(firstLine);
        n2 = numel(secondLine);
        len = max(n1, n2);
        
        lim = [repmat('*', 1, len+4), '\n'];
        
        if i > 1
            fprintf('\n')
        end
        fprintf(lim);
        printstr = ['*%-', num2str(len+1), 's *\n'];
        fprintf(printstr, firstLine);
        fprintf(printstr, secondLine);
        fprintf(lim);
        [state0, model, schedule, restartStep, restartOffset, msg] = getRestart(problem, opt);
        if isnan(restartStep)
            if opt.continueOnError
                fprintf(msg);
            else
                error(msg);
            end
        end
        timer = tic();
        if isfinite(restartStep)
            mods = mrstModule();
            try
                mrstModule('add', problem.Modules{:});
                simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls,...
                    'restartStep', restartStep,...
                    'OutputHandler', state_handler, ...
                    'WellOutputHandler', wellSol_handler, ...
                    'ReportHandler', report_handler, ...
                    'outputOffset', restartOffset, ...
                    problem.SimulatorSetup.ExtraArguments{:}, ...
                    'OutputMinisteps', doMinisteps);
            catch ex
                mrstModule('reset', mods{:});
                msg = ex.message;
                ok = false;
                fprintf('!!! Simulation resulted in fatal error !!!\n Exception thrown: %s\n', msg);
                if ~opt.continueOnError
                    rethrow(ex);
                end
            end
            if ok
                mrstModule('reset', problem.Modules{:});
            end
        end
        
        status(i).Message = msg;
        status(i).WallTime = toc(timer);
        status(i).RestartIndex = restartStep;
    end
    if opt.plot
        plotPackedProblem(problem);
    end
end

function [state0, model, schedule, restart, outoffset, msg] = getRestart(problem, opt)
    schedule = problem.SimulatorSetup.schedule;
    model = problem.SimulatorSetup.model;
    state0 = problem.SimulatorSetup.state0;
    T = cumsum(schedule.step.val);
    
    state_handler = problem.OutputHandlers.states;
    doMinisteps = problem.SimulatorSetup.OutputMinisteps;
    nstep = numel(schedule.step.val);
    ndata = state_handler.numelData();
    msg = '';
    outoffset = [];
    if opt.checkTooMany && ~doMinisteps
        % This is a serious error!
        assert(ndata <= nstep, 'Too much data exists for %s! Problem may have been redefined.', problem.Name);
    end
    restart = opt.restartStep;
    if ndata == 0 || restart <= 1
        endstate = state0;
        state0.time = 0;
    else
        endstate = state_handler{ndata};
    end
    autoRestart = isnan(restart);
    if autoRestart
        % No restart specified. We need to figure out where to pick up the
        % simulation, or if it is already done, determined from whatever
        % data is available.
        if ndata == 0
            % We start from the beginning
            restart = 1;
        elseif (ndata >= nstep && ~doMinisteps) || (endstate.time >= T(end) && doMinisteps)
            % The simulation is already done
            restart = inf;
        else
            % We have a partial simulation. Find last complete step, and
            % then grab that state as the restart.
            if doMinisteps
                % Find last state, then find the corresponding control step
                last = find(T <= endstate.time, 1, 'last');
                if isempty(last)
                    restart = 1;
                else
                    restart = last + 1;
                    [state0, outoffset] = state_handler.getByTime(T(last));
                    outoffset = outoffset + 1;
                end
            else
                restart = ndata + 1;
                state0 = endstate;
            end
        end
    elseif restart > 1
        % Requested specific restart step. We only need to do something if
        % restart > 1. Otherwise, state0 + restart = 1 is fine.
        T_prev = T(restart-1);
        if endstate.time < T_prev
            % This is bad!
            msg = sprintf('Restart step %d was specified, but step %d is not complete! Aborting.', restart, restart-1);
            restart = nan;
        elseif doMinisteps
            [state0, outoffset] = state_handler.getByTime(T_prev);
            outoffset = outoffset + 1;
        else
            state0 = state_handler{restart-1};
        end
    end
    % Might happen for floating point comparisons of ministep time?
    if isempty(state0)
        restart = nan;
        msg = 'Internal error. Restart state not found.';
    end
    % Print out some messages
    if isinf(restart)
        fprintf('-> Complete output found, nothing to do here.\n');
    elseif restart == 1
        if autoRestart
            fprintf('-> No output found. Starting from first step...\n');
        else
            fprintf('-> Starting from first step as requested...\n');
        end
    else
        if autoRestart
            fprintf('-> Partial output found. Starting from step %d of %d...\n', restart, nstep);
        else
            fprintf('-> Partial output found. Starting from requested step %d of %d...\n', restart, nstep);
        end
    end
end
