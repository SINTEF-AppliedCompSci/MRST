function [wellSols, states] = runScheduleRefactor(initState, model, schedule, varargin)

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

    wellSols = cell(nSteps, 1);
    states   = cell(nSteps, 1);


    getWell = @(index) schedule.control(schedule.step.control(index)).W;
    state = initState;
    if ~isfield(state, 'wellSol')
        state.wellSol = initWellSolLocal(getWell(1), state);
    end

    for i = 1:nSteps
        fprintf('Solving timestep %d of %d at %s\n', i, nSteps, formatTimeRange(tm(i)));
        W = getWell(i);
        timer = tic();
        [state, status] = solver.solveTimestep(state, dt(i), model, 'Wells', W);
        t = toc(timer);
        dispif(vb, 'Completed %d iterations in %2.2f seconds (%2.2fs per iteration)\n', ...
                    status.iterations, t, t/status.iterations);
        states{i} = state;
        wellSols{i} = state.wellSol;
    end
end
%--------------------------------------------------------------------------

