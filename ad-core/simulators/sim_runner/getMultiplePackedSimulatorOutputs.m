function [all_ws, all_states, all_reports, names, T, grids, act] = getMultiplePackedSimulatorOutputs(problems, varargin)
    opt = struct('minSteps', 0);
    [opt, arg] = merge_options(opt, varargin{:});
    np = numel(problems);
    counts = zeros(np, 1);
    [T, all_states, all_ws, all_reports, names, grids] = deal(cell(np, 1));
    for i = 1:np
        problem = problems{i};
        [all_ws{i}, all_states{i}, all_reports{i}] = getPackedSimulatorOutput(problem, arg{:});
        names{i} = problem.Name;
        counts(i) = numel(all_ws{i});
        T{i} = cumsum(problem.SimulatorSetup.schedule.step.val(1:counts(i)));
        grids{i} = problem.SimulatorSetup.model.G;
    end
    
    act = counts >= opt.minSteps;
    
    all_ws = all_ws(act);
    all_states = all_states(act);
    all_reports = all_reports(act);
    names = names(act);
    T = T(act);
    grids = grids(act);
end