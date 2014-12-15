function solver = getNonLinearSolver(model, varargin)
% Set up reasonable defaults for the nonlinear solver for a field
% simulation with significant size and complexity

opt = struct('DynamicTimesteps', true, ...
             'useCPR',           true);

[opt, varg] = merge_options(opt, varargin{:});

[linsolve, timestepper] = deal([]);
if opt.useCPR && isa(model, 'ReservoirModel')
    if ~isempty(mrstPath('agmg'))
        mrstModule add agmg
        pSolver = AGMGSolverAD();
    else
        pSolver = BackslashSolverAD();
    end
    linsolve = CPRSolverAD('ellipticSolver', pSolver);
end

if opt.DynamicTimesteps
    timestepper = ...
    IterationCountTimeStepSelector('firstRampupStepRelative', 0.1, ...
                                   'firstRampupStep',         1*day);
end
solver = NonLinearSolver('timeStepSelector', timestepper, ...
                         'LinearSolver', linsolve, varg{:});
end