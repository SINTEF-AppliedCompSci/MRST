function nls = getDefaultFlashNonLinearSolver()
% Get default nonlinear solver for flash problems.
    nls = NonLinearSolver();
    nls.verbose = -1;
    nls.maxIterations = 100;
    nls.maxTimestepCuts = 0;
    nls.useRelaxation = true;
    nls.continueOnFailure = true;
    nls.errorOnFailure = false;
end