function [states, failure] = runSimulationProblem(model, state0, schedule, varargin)

    opt = struct('useCPR',  false, ...
                 'useAGMG', false, ...
                 'stepcount',   inf, ...
                 'cprTol', 1e-2 ...
                );
   
    opt = merge_options(opt, varargin{:});
    
    mrstModule add ad-fi ad-refactor
    
    if isfinite(opt.stepcount)
        schedule.step.val = schedule.step.val(1:opt.stepcount);
        schedule.step.control = schedule.step.control(1:opt.stepcount);
    end
    
    if opt.useCPR
        if opt.useAGMG
            mrstModule add agmg
            ellipSolver = AGMGSolverAD();
        else
            ellipSolver = BackslashSolverAD();
        end
        linsolve = CPRSolverAD('ellipticSolver', ellipSolver, ...
                               'relativeTolerance', opt.cprTol);
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