function problem = initEclipsePackedProblemAD(deck, varargin)
    opt = struct('name', []);
    [opt, extra] = merge_options(opt, varargin{:});
    
    [state0, model, schedule, nls] = initEclipseProblemAD(deck, extra{:});
    if isempty(opt.name)
        name = model.inputdata.RUNSPEC.TITLE;
    else
        name = opt.name;
    end
    problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);
end