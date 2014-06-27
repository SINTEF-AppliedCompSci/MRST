function [wellSols, states, schedulereport] = runScheduleRefactor(initState, model, schedule, varargin)

    opt = struct('Verbose', mrstVerbose,...
                 'linearSolver', []);

    opt = merge_options(opt, varargin{:});

    vb = opt.Verbose;
    %--------------------------------------------------------------------------

    dt = schedule.step.val;
    tm = cumsum(dt);
    dispif(vb, '*****************************************************************\n')
    dispif(vb, '********** Starting simulation: %5.0f steps, %5.0f days *********\n', numel(dt), tm(end)/day)
    dispif(vb, '*****************************************************************\n')

    solver = nonlinearSolver('linearSolver', opt.linearSolver);

    nSteps = numel(dt);

    [wellSols, states, reports] = deal(cell(nSteps, 1));
    wantStates = nargout > 1;
    wantReport = nargout > 2;

    getWell = @(index) schedule.control(schedule.step.control(index)).W;
    state = initState;
    if ~isfield(state, 'wellSol')
        state.wellSol = initWellSolLocal(getWell(1), state);
    end
    
    
    simtime = zeros(nSteps, 1);
    for i = 1:nSteps
        fprintf('Solving timestep %d of %d at %s\n', i, nSteps, formatTimeRange(tm(i)));
        W = getWell(i);
        timer = tic();
        
        [state, report] = solver.solveTimestep(state, dt(i), model, 'Wells', W);
        t = toc(timer);
        dispif(vb, 'Completed %d iterations in %2.2f seconds (%2.2fs per iteration)\n', ...
                    report.Iterations, t, t/report.Iterations);
        
        wellSols{i} = state.wellSol;
        
        if wantStates
            states{i} = state;
        end
        
        if wantReport
            reports{i} = report;
        end
    end
    
    if wantReport
        schedulereport = struct();
        schedulereport.ControlstepReports = reports;
        schedulereport.ReservoirTime = cumsum(schedule.step.val);
        schedulereport.Converged  = cellfun(@(x) x.Converged, reports);
        schedulereport.Iterations = cellfun(@(x) x.Iterations, reports);
        schedulereport.SimulationTime = simtime;
    end
end
%--------------------------------------------------------------------------

