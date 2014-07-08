function [wellSols, states, schedulereport] = simulateScheduleAD(initState, model, schedule, varargin)

    opt = struct('Verbose', mrstVerbose,...
                 'OutputMinisteps', false, ...
                 'NonLinearSolver', [], ...
                 'OutputHandler',   [], ...
                 'LinearSolver', []);

    opt = merge_options(opt, varargin{:});

    vb = opt.Verbose;
    %--------------------------------------------------------------------------

    dt = schedule.step.val;
    tm = cumsum(dt);
    dispif(vb, '*****************************************************************\n')
    dispif(vb, '********** Starting simulation: %5.0f steps, %5.0f days *********\n', numel(dt), tm(end)/day)
    dispif(vb, '*****************************************************************\n')
    
    solver = opt.NonLinearSolver;
    if isempty(solver)
        solver = NonLinearSolver('linearSolver', opt.LinearSolver);
    elseif ~isempty(opt.LinearSolver)
        % We got a nonlinear solver, but we still want to override the
        % actual linear solver passed to the higher level schedule function
        % we're currently in
        solver.LinearSolver = opt.LinearSolver;
    end
    nSteps = numel(dt);

    [wellSols, states, reports] = deal(cell(nSteps, 1));
    wantStates = nargout > 1;
    wantReport = nargout > 2;

    getWell = @(index) schedule.control(schedule.step.control(index)).W;
    state = initState;
    if ~isfield(state, 'wellSol')
        state.wellSol = initWellSolLocal(getWell(1), state);
    end
    
    failure = false;
    simtime = zeros(nSteps, 1);
    prevControl = nan;
    for i = 1:nSteps
        fprintf('Solving timestep %d of %d at %s\n', i, nSteps, formatTimeRange(tm(i)));
        currControl = schedule.step.control(i);
        if prevControl ~= currControl 
            W = schedule.control(currControl).W;
            forces = model.getDrivingForces(schedule.control(currControl));
            prevControl = currControl;
        end

        timer = tic();
        
        
        state0 = state;
        state0.wellSol = initWellSolLocal(W, state);
        
        if opt.OutputMinisteps
            [state, report, ministeps] = solver.solveTimestep(state0, dt(i), model, ...
                                            forces{:}, 'controlId', currControl);
        else
            [state, report] = solver.solveTimestep(state0, dt(i), model,...
                                            forces{:}, 'controlId', currControl);
        end
        t = toc(timer);
        
        if ~report.Converged
            warning('Nonlinear solver aborted, returning incomplete results!');
            failure = true;
            break;
        end
        dispif(vb, 'Completed %d iterations in %2.2f seconds (%2.2fs per iteration)\n', ...
                    report.Iterations, t, t/report.Iterations);
        

        W = updateSwitchedControls(state.wellSol, W);
        
        
        % Handle massaging of output to correct expectation
        if opt.OutputMinisteps
            % We have potentially several ministeps desired as output
            nmini = numel(ministeps);
            ise = find(cellfun(@isempty, states), 1, 'first');
            if isempty(ise)
                ise = numel(states) + 1;
            end
            ind = ise:(ise + nmini - 1);
            states_step = ministeps;
        else
            % We just want the control step
            ind = i;
            states_step = {state};
        end
        wellSols_step = cellfun(@(x) x.wellSol, states_step, 'UniformOutput', false);
        
        wellSols(ind) = wellSols_step;
        
        if ~isempty(opt.OutputHandler)
            opt.OutputHandler{ind} = states_step;
        end
        
        if wantStates
            states(ind) = states_step;
        end
        
        if wantReport
            reports{i} = report;
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
end
%--------------------------------------------------------------------------

