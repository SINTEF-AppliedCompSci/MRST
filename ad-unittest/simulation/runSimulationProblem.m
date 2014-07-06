function [states, failure] = runSimulationProblem(model, state0, schedule, varargin)

    opt = struct('useCPR',  false, ...
                 'useAGMG', false, ...
                 'cprTol', 1e-2 ...
                );
   
    opt = merge_options(opt, varargin{:});
    
    mrstModule add ad-fi ad-refactor
    
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
                                'maxSubsteps', 2,...
                                'errorOnFailure', false);
    
    [wellSols, states, report] = simulateScheduleAD(state0, model, schedule,...
                                'NonLinearSolver', nonlinear);
    failure = report.Failure;
end