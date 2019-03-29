function [states, failure] = runSimulationProblem(model, state0, schedule, varargin)

    opt = struct('useCPR',  false, ...
                 'useAGMG', false, ...
                 'selectLinearSolver', false, ...
                 'stepcount',   inf, ...
                 'cprTol', 1e-3 ...
                );
   
    opt = merge_options(opt, varargin{:});
    
    if isfinite(opt.stepcount)
        schedule.step.val = schedule.step.val(1:opt.stepcount);
        schedule.step.control = schedule.step.control(1:opt.stepcount);
    end
    
    if opt.useCPR
        if opt.useAGMG
            mrstModule add agmg
            ellipSolver = AGMGSolverAD('tolerance', 1e-2');
        else
            ellipSolver = BackslashSolverAD();
        end
        linsolve = CPRSolverAD('ellipticSolver', ellipSolver, ...
                               'relativeTolerance', opt.cprTol);
    elseif opt.selectLinearSolver
        linsolve = selectLinearSolverAD(model);
    else
        linsolve = BackslashSolverAD();
    end
    % Set max substeps low because it is a failure if these tests do not
    % run through the entire schedule without cutting timesteps more than
    % once
    nonlinear = NonLinearSolver('LinearSolver', linsolve, ...
                                'maxTimestepCuts', 2,...
                                'errorOnFailure', false);
    % Ignore flux output for tests
    model.outputFluxes = false;
    
    [wellSols, states, report] = simulateScheduleAD(state0, model, schedule,...
                                'NonLinearSolver', nonlinear);
    failure = report.Failure;
end