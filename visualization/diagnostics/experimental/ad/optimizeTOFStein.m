function [D_best W_best history] = optimizeTOFStein(G, W, mob, pv, T, N, s, state0, minvals, objective, varargin)
opt = struct('maxiter', 25, ...
             'alpha', 1, ...
             'targets', [], ...
             'deltatol', 1e-6, ...
             'plotProgress', false, ...
             'verbose', mrstVerbose, ...
             'linsolve', @mldivide, ...
             'alphamod', 10);

opt = merge_options(opt, varargin{:});

[state, D, grad] = solveStationaryPressureStein(G, state0, s, W, mob, pv, T, N, 'objective', objective, 'linsolve', opt.linsolve);
if isempty(opt.targets)
    targets = D.inj;
else
    targets = opt.targets;
end

assert(numel(targets) == numel(minvals) || numel(minvals) == 1);
targets = reshape(targets, [], 1);
minvals = reshape(minvals, [], 1);
assert(~any(ismember(targets, D.inj)) || ~any(ismember(targets, D.prod)),...
    'Optimizing both injectors and producers at the same time is not possible')

W_previous = W;
W_best = W;
D_best = D;
D_all = D;
grad_all = grad;
W_vals = [W(targets).val];

iterhist = 0;

% Baseline objective function
obj0 = grad.objective.val;
objval = obj0;

% Scale initial alpha value
alpha_0 = abs(opt.alpha*grad.objective.val/max(abs(grad.well)));
alpha = alpha_0;

dispif(opt.verbose, 'Obj: %2.6g at initial controls \n', objval);
outeriter = 0;
while 1
    improved = false;
    inneriter = 0;
    outeriter = outeriter + 1;
    while ~improved
        values = vertcat(W(targets).val);

        % Use signs here to ensure that the step selector does not have to care
        % about them...
        wsign = 1 - 2*ismember(targets, D.prod);

        g = wsign.*grad.well(targets);
        g = wsign.*projectGradient(g, wsign.*values, alpha, wsign.*minvals, []);
        values = update(values, alpha, g, []);
        for i = 1:numel(targets)
            W(targets(i)).val = values(i);
        end
        W_vals = [W_vals; values .']; %#ok

        if opt.plotProgress
            plot(W_vals)
            drawnow
        end


        [state, D, gradnew] = solveStationaryPressureStein(G, state0, s, W, mob, pv, T, N, 'objective', objective);
        newobj = gradnew.objective.val;
        %gradnew.well


        delta = (newobj-objval)/objval;
        dispif(opt.verbose, 'Obj: %2.6g (delta: %1.1g, alpha: %1.1g) at iteration %d\n', newobj, delta, alpha, inneriter);
        inneriter = inneriter + 1;



        if newobj < objval && abs(delta) > opt.deltatol
            dispif(opt.verbose, '================ VALUE IMPROVED ==================\n')
            alpha = alpha*opt.alphamod;
            improved = true;
        end

        if improved || inneriter > opt.maxiter  || abs(delta) < opt.deltatol || any(~isfinite(g))
            break;
        end

        alpha = alpha/opt.alphamod;
        W = W_previous;

    end
    iterhist = [iterhist inneriter]; %#ok
    if ~improved
        break
    end

    objval = newobj;
    W_best = W;
    D_best = D;
    D_all = [D_all D]; %#ok
    grad_all = [grad_all; gradnew]; %#ok
    W_previous = W;
    grad = gradnew;
end

history.values = W_vals;
history.D = D_all;
history.gradient = grad_all;
history.iterations = iterhist + 1;
dispif(opt.verbose, 'Optimization loop DONE, improvement %2.2g (relative improvement: %2.2g)\n', obj0 - newobj, abs(obj0-newobj)/obj0);
end


function [values step] = update(values, alpha, grad, objective)
    if 1
        step = alpha*grad;
        values = values - step;
    else
        step = objective./grad;
        values = values - step;
    end
end


function grad = projectGradient(grad, values, alpha, lowlimit, obj)
    iter = 0;
    grad = grad - mean(grad);
    while 1

        [values_new step] = update(values, alpha, grad, obj);
        overshoot = abs(min(values_new - lowlimit, 0));
        iter = iter + 1;

        adjustment = 1 - min(overshoot./abs(step), 1);
        grad = grad.*adjustment;

        grad = grad - mean(grad);

        if norm(overshoot, inf) < 1e-8 && abs(sum(grad))/max(abs(grad)) < 1e-8
            break
        end
        if iter > 10000
            warning('Unable to project gradient to constraints in 1000 iterations...')
            break
        end
    end
end
