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
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
        schedule = problem.SimulatorSetup.schedule;
        state0 = problem.SimulatorSetup.state0;
        model = problem.SimulatorSetup.model;
        nls = problem.SimulatorSetup.NonLinearSolver;
        
        state_handler = problem.OutputHandlers.states;
        wellSol_handler = problem.OutputHandlers.wellSols;
        report_handler = problem.OutputHandlers.reports;

        nstep = numel(schedule.step.val);
        ndata = state_handler.numelData();
        ok = true;
        doSim = true;
        msg = '';
        if opt.checkTooMany
            % This is a serious error!
            assert(ndata <= nstep, 'Too much data exists for %s! Problem may have been redefined.', problems{i}.Name);
        end
        
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
        if ~isnan(opt.restartStep)
            restartStep = opt.restartStep;
            if ndata < restartStep - 1
                msg = sprintf('Restart step %d was specified, but step %d is not complete! Aborting.', restartStep, restartStep-1);
                doSim = false;
                if ~opt.continueOnError
                    error(msg); %#ok
                end
            end
            if restartStep > 1
                state0 = state_handler{restartStep-1};
            end
        elseif ndata == nstep
            fprintf('-> Complete output found, nothing to do here.\n');
            % Already run!
            doSim = false;
            restartStep = nan;
        elseif ndata == 0
            fprintf('-> No output found, starting from first step...\n');
            restartStep = 1;
        else
            fprintf('-> Partial output found, starting from step %d of %d...\n', ndata+1, nstep);
            state0 = state_handler{ndata};
            restartStep = ndata + 1;
        end

        timer = tic();
        if doSim
            mods = mrstModule();
            try
                mrstModule('add', problem.Modules{:});
                simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls,...
                                                            'restartStep', restartStep,...
                                                            'OutputHandler', state_handler, ...
                                                            'WellOutputHandler', wellSol_handler, ...
                                                            'ReportHandler', report_handler, ...
                                                            problem.SimulatorSetup.ExtraArguments{:});
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