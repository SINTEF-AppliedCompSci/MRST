function [ws, states] = simulatePackedProblemJutul(problem, varargin)
    s = problem.SimulatorSetup;
    state0 = s.state0;
    schedule = s.schedule;
    model = s.model;
    name = [problem.BaseName, '_', problem.Name];
    name = name(~isspace(name));
    [ws, states] = simulateScheduleJutul(state0, model, schedule, 'name', name, varargin{:});
end
