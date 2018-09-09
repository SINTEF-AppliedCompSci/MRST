function [all_ws, all_states, all_reports, names, T, grids] = getMultiplePackedSimulatorOutputs(problems, varargin)
    np = numel(problems);
    [T, all_states, all_ws, all_reports, names, grids] = deal(cell(np, 1));
    for i = 1:np
        problem = problems{i};
        [all_ws{i}, all_states{i}, all_reports{i}] = getPackedSimulatorOutput(problem, varargin{:});
        names{i} = problem.Name;
        n = numel(all_ws{i});
        T{i} = cumsum(problem.SimulatorSetup.schedule.step.val(1:n));
        grids{i} = problem.SimulatorSetup.model.G;
    end
end