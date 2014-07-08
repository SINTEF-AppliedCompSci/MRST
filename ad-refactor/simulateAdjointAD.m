function gradients = simulateAdjointAD(state0, states, model, schedule, getObjective, varargin)
    opt = struct('ControlVariables', 'well', ...
                 'LinearSolver', []);
    opt = merge_options(opt, varargin{:});
    
    getState = @(i) getStateFromCell(states, state0, i);
    
    if isempty(opt.LinearSolver)
        linsolve = BackslashSolverAD();
    else
        linsolve = opt.LinearSolver;
    end
    
    if iscell(opt.ControlVariables)
        % Path presently not supported by merge_options limitations
        ncv = numel(opt.ControlVariables);
    else
        ncv = 1;
    end
    nstep = numel(schedule.step.val);
    grad = [];
    gradstep = cell(nstep, ncv);
    nt = nstep;
    for step = nt:-1:1
        [dg, grad, report] = model.solveAdjoint(linsolve, getState, ...
                                                getObjective, schedule, grad, step);
        gradstep(step, :) = getRequestedGradients(dg, report, opt.ControlVariables);
    end
    
    
    % Sum up to the control steps
    nc = numel(schedule.control);
    gradients = cell(nc, ncv);
    for k = 1:nc
        ck = schedule.step.control == k;
        for j = 1:ncv
            tmp = gradstep(ck, j);
            gradients{k, j} = sum(horzcat(tmp{:}), 2);
        end
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

function g = getRequestedGradients(dg, report, wantGradFor)
    if ischar(wantGradFor)
        g = {vertcat(dg{strcmpi(report.Types, wantGradFor)})};
    else
        ng = numel(wantGradFor);
        g = cell(1, ng);
        for i = 1:ng
            n = wantGradFor{i};
            g{i} = vertcat(dg{strcmpi(report.Types, n)});
        end
    end
end