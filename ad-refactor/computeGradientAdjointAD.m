function gradients = computeGradientAdjointAD(state0, states, model, schedule, getObjective, varargin)
    opt = struct('ControlVariables', 'well', ...
                 'Scaling',          [], ...
                 'LinearSolver',     []);
    opt = merge_options(opt, varargin{:});
    
    getState = @(i) getStateFromInput(schedule, states, state0, i);
    
    if isempty(opt.LinearSolver)
        linsolve = BackslashSolverAD();
    else
        linsolve = opt.LinearSolver;
    end
    
    if isempty(opt.Scaling)
        opt.Scaling = struct('rate', 1, 'pressure', 1);
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
                                         getObjective, schedule, grad, step, opt.Scaling);
        gradstep(step, :) = getRequestedGradients(dg, report, opt.ControlVariables);
    end
    
    
    % Sum up to the control steps
    nc = numel(schedule.control);
    gradients = cell(ncv, nc);
    for k = 1:nc
        ck = schedule.step.control == k;
        for j = 1:ncv
            tmp = gradstep(ck, j);
            gradients{j, k} = full(sum(horzcat(tmp{:}), 2));
        end
    end
end

function state = getStateFromInput(schedule, states, state0, i)
    if i == 0
        state = state0;
    elseif i > numel(schedule.step.val)
        state = [];
    else
        state = states{i};
    end
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