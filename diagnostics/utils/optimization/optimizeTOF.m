function [D_best, W_best, history] = optimizeTOF(G, W, fluid, pv, T, op, state, minRates, objective, varargin)
%Optimize well rates based on diagnostics-based objective function
%
% SYNOPSIS:
%  [D_best, W_best, history] = optimizeTOF(G, W, fluid, pv, T, s, state, minRate, objective)
%  [D_best, W_best, history] = optimizeTOF(G, W, fluid, pv, T, s, state, minRate, objective, 'targets', [1, 3, 9])
%
% DESCRIPTION:
%   Optimizes well rates. The optimization process is based on steepest
%   descent given a gradient, and uses the following constraints to perform
%   the optimization:
%
%   - Sum of rates of the targets are constant, i.e. the total
%   injected/produced volume of the targets remains constant throughout the
%   optimization.
%
%   - All target wells must have a minimum rate which ensures that they do
%   not change from injectors to producers during the optimization process.
%   This can be specified for all wells simultanously, or well specific
%   limits can be provided.
%
%   The progress of the optimization can be visualized during the process
%   (see optional parameters below)
%
%
% REQUIRED PARAMETERS:
%   G     - Valid grid structure.
%
%   W     - Well configuration suitable for AD-solvers.
%
%   fluid - AD-fluid, for instance defined by initSimpleADIFluid. This is
%           used to evalute relperm / viscosity when calculating fluxes.
%
%   pv    - Pore volume. Typicall output from poreVolume(G, rock).
%
%   op    - Operators from setupOperatorsTPFA
%
%   state - Reservoir state as defined by initResSol.
%
%   minRates - Either a single value for the minimum well rate, or a value
%              per target well (see 'targets' keyword)
%
%   objective - Objective function handle. Should have the interface
%               @(state, D) where D is a diagnostics object (see
%               computeTOFandTracer for the expected fields) and return a
%               single numerical value along with its Jacobian. See
%               getObjectiveDiagnostics for examples.
%
%
% OPTIONAL PARAMETERS:
%
%  maxiter - Maximum number of attempts to find an improvement. If maxiter
%            successive evaluations fail to improve upon the last optimum,
%            the function returns.
%
%  alpha   - Adjustment for scaling the initial relaxation factor in the
%            optimization algorithm. This can improve convergence if the
%            problem if chosen wisely, but the default is fairly robust
%            since it is based on the magnitude of the initial objective as
%            well as the gradients of the well at i = 0.
%
%  targets - Index into W argument of wells that are to be optimized. They
%            *must* all be initially the same type (i.e. producers /
%            injectors) and they should be rate controlled. This defaults
%            to all injectors with rate controls.
%
%  deltatol - Termination criterion. If the relative change between
%            iterations falls below this number, the function returns. For
%            instance, if deltatol is 0.05, if the relative change between
%            evaluations is less than 5% in absolute value, the iteration
%            completes.
%
%  plotProgress - Plots the well rates and function evaluations during the
%                 well optimization process.
%
%  verbose   - Controls the amount of output produced by the routine.
%
%  linsolve  - The linear solver used to solve the elliptic pressure
%              systems during the optimization with a interface on the form
%              @(A, b) for the linear system Ax=b. Fast multigrid or
%              multiscale methods will significantly improve upon the
%              solution speed for larger problems.
%
%  alphamod  - The modifier used to adjust alpha during the process. If a
%              step is unsuccessful,
%                     alpha <- alpha/alphamod
%              and the step is retried. If a step is successful, the next
%              initial alpha is set to
%                     alpha <- alpha*alphamod.
%
% RETURNS:
%   D_best   - The flow diagnostics object corresponding to the best
%              evaluated well configuration.
%
%   W_best   - The best well configuration
%
%   history  - A struct containing the optimization history in the fields:
%                   D - Diagnostic objects for all steps that reduced the
%                   objective function.
%
%                   values - Values per target well. One row per element in
%                            the D array.
%
%                   gradient - One gradient per diagnostics object. The
%                              objective function value for step i can be
%                              found at history.gradient(1).objective.val.
%
%                   iterations - Iteration per step required to improve
%                               upon the previous minimum.


%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


opt = struct('maxiter', 25, ...
             'alpha', 1, ...
             'targets', [], ...
             'deltatol', 1e-3, ...
             'plotProgress', false, ...
             'verbose', mrstVerbose, ...
             'linsolve', @mldivide, ...
             'msbasis', [], ...
             'autoscalealpha', true, ...
             'alphamod', 10);

opt = merge_options(opt, varargin{:});

[st, D, grad] = solveStationaryPressure(G, state, op, W, fluid, pv, T, 'objective', objective,...
                    'linsolve', opt.linsolve, 'msbasis', opt.msbasis);
if isempty(opt.targets)
    tmp = false(numel(W), 1);
    tmp(D.inj) = true;
    % Pick all rate controlled injectors as default
    targets = find(tmp & arrayfun(@(x) strcmpi(x.type, 'rate'), W));
else
    targets = opt.targets;
end

assert(numel(targets) == numel(minRates) || numel(minRates) == 1);
targets = reshape(targets, [], 1);
minRates = reshape(minRates, [], 1);
assert(~any(ismember(targets, D.inj)) || ~any(ismember(targets, D.prod)),...
    'Optimizing both injectors and producers at the same time is not possible')
assert(all(arrayfun(@(x) strcmpi(x.type, 'rate'), W(targets))), ...
    'optimizeTOF is only tested for rate based wells!');
W_previous = W;
W_best = W;
D_best = D;
D_all = D;
grad_all = grad;
W_vals = [W(targets).val];
improvements = [];
iterhist = 0;

% Baseline objective function
obj0 = grad.objective.val;
objval = obj0;
objectivevalues = [];

% Scale initial alpha value
if opt.autoscalealpha
    alpha = abs(opt.alpha*grad.objective.val/max(abs(grad.well)));
else
    alpha = opt.alpha;
end

dispif(opt.verbose, 'Obj: %2.6g at initial controls \n', objval);
outeriter = 0;
while true
    improved = false;
    inneriter = 0;
    outeriter = outeriter + 1;
    while ~improved
        values = vertcat(W(targets).val);

        % Use signs here to ensure that the step selector does not have to care
        % about them...
        wsign = 1 - 2*ismember(targets, D.prod);

        g = wsign.*grad.well(targets);
        g = wsign.*projectGradient(g, wsign.*values, alpha, wsign.*minRates);
        values = update(values, alpha, g);
        for i = 1:numel(targets)
            W(targets(i)).val = values(i);
        end
        W_vals = [W_vals; values .']; %#ok

        % Calculate time of flight and gradients
        [st, D, gradnew] = solveStationaryPressure(G, state, ...
           op, W, fluid, pv, T, 'objective', objective, 'linsolve', opt.linsolve);
        newobj = gradnew.objective.val;


        delta = (newobj-objval)/objval;
        dispif(opt.verbose, 'Obj: %2.6g (delta: %1.1g, alpha: %1.1g) at iteration %d\n', newobj, delta, alpha, inneriter);
        inneriter = inneriter + 1;

        if newobj < objval && abs(delta) > opt.deltatol
            dispif(opt.verbose, '================ VALUE IMPROVED ==================\n')
            alpha = alpha*opt.alphamod;
            improved = true;
        end
        improvements = [improvements; improved]; %#ok
        objectivevalues = [objectivevalues; newobj];%#ok

        if opt.plotProgress
            % Draw the well rates
            clf;
            subplot(2,1,1)
            plotWellRates(W(targets), W_vals, improvements > 0)
            title('Well rates')
            xlabel('Optimization step')
            ylabel('Normalized injection volume')

            % And draw the objective values
            subplot(2,1,2);
            hold on
            stairs(objectivevalues, '--', 'linewidth', 2)
            tmp = find(improvements > 0);
            stairs(tmp, objectivevalues(tmp), 'd', 'markerfacecolor', 'g', 'markersize', 12)
            title('Objective function value');
            grid on
            grid minor
            axis tight
            ylim([min(0, min(objectivevalues)), max(objectivevalues)])
            drawnow
        end

        if improved || inneriter > opt.maxiter  || abs(delta) < opt.deltatol || any(~isfinite(g))
            % We either have an improved solution or we are giving up
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


function [values, step] = update(values, alpha, grad)
    step = alpha*grad;
    values = values - step;
end


function grad = projectGradient(grad, values, alpha, lowlimit)

    iter = 0;
    % Ensure total fluid volume
    grad = grad - mean(grad);
    while 1
        [values_new, step] = update(values, alpha, grad);
        overshoot = abs(min(values_new - lowlimit, 0));
        iter = iter + 1;

        adjustment = 1 - min(overshoot./abs(step), 1);
        grad = grad.*adjustment;

        grad = grad - mean(grad);

        if norm(overshoot, inf) < 1e-8 && abs(sum(grad))/max(abs(grad)) < 1e-8
            break
        end
        if iter > 10000
            warning(['Unable to project gradient to constraints in 10000'...
                    ' iterations. Total volume injection may be violated.']);
            break
        end
    end
end
