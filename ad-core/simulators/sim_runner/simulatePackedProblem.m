function [ok, status] = simulatePackedProblem(problems, varargin)
    opt = struct('checkTooMany', true, ...
                 'continueOnError', true);
    opt = merge_options(opt, varargin{:});
    if isstruct(problems)
        problems = {problems};
    end
    
    assert(iscell(problems));
    
    np = numel(problems);
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
                problem.Description, problem.Name);
        
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
        if ndata == nstep
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
end