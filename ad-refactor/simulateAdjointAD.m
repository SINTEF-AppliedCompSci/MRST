function grad = simulateAdjointAD(state0, states, model, schedule, getObjective, varargin)
    
    getState = @(i) getStateFromCell(states, state0, i);
    
    linsolve = BackslashSolverAD();
    grad = [];
    
    nt = numel(schedule.step.val);
    for step = nt:-1:1
        [grad, report] = model.solveAdjoint(linsolve, getState, ...
                                                getObjective, schedule, grad, step);
    end
end

function state = getStateFromCell(states, state0, i)
    if i == 0
        state = state0;
    elseif i > numel(states)
        state = [];
    else
        state = states{i};
    end
end

function state = getStateFromFile(i)
    assert(0, 'Someone should implement me!');
end