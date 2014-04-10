function [wellSols, states] = runScheduleRefactor(initState, model, schedule, varargin)

    opt = struct('Verbose', mrstVerbose);

    opt = merge_options(opt, varargin{:});

    vb = opt.Verbose;
    %--------------------------------------------------------------------------

    dt = schedule.step.val;
    tm = cumsum(dt);
    dispif(vb, '*****************************************************************\n')
    dispif(vb, '********** Starting simulation: %5.0f steps, %5.0f days *********\n', numel(dt), tm(end)/day)
    dispif(vb, '*****************************************************************\n')

    solver = nonlinearSolver();

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
        [state, status] = solver.solveTimestep(state, dt(i), model, 'Wells', W);
        states{i} = state;
        wellSols{i} = state.wellSol;
    end
end
%--------------------------------------------------------------------------

