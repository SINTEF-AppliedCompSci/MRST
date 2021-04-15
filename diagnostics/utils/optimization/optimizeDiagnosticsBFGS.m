function [states, Ds, info] = optimizeDiagnosticsBFGS(G, W, fluid, pv, T, s, state, boxLims, objective, varargin)
% Optimize well rates/bhps to minimize/maximize objective 
%
% SYNOPSIS:
%  [W_best, D_init, D_best] = optimizeDiagnosticsBFGS(G, W, fluid, pv, T, s, state, minRate, objective)
%
% DESCRIPTION:
%   Optimizes well rates/bhps as defined in W. The optimization process is based 
%   on BFGS where the search direction is projected onto the feasible domain
%   as given by boxLims (upper/lower value for each control). 
%   
%   - boxLims must be supplied for all target wells.
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
%   s     - AD system. See initADISystem.
%
%   state - Reservoir state as defined by initResSol.
%
%   boxLims  - Lower and upper limit for each target well 
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

opt = {'targets', (1:numel(W)).', ...
       'mbasis' ,             [], ...
       'computeBasis',     false, ...
       'objectiveScaling',    [], ...
       'reportTimings',     true, ...
       'linsolve',     @mldivide, ...
       'minimize',          true};

optBFGS = {'gradTol',             1e-3, ...
           'objChangeTol',        1e-4, ...
           'maxIt',               100,   ...
           'lineSearchMaxIt',     5,   ...
           'stepInit',            -1,   ...
           'wolfe1',              1e-3, ...
           'wolfe2',              0.9,  ...
           'safeguardFac',        0.0001, ...  
           'stepIncreaseTol',     10,    ...
           'useBFGS',             true, ...
           'linEq',                 []};
       
opt = struct(opt{:}, optBFGS{:});
opt = merge_options(opt, varargin{:});

% setup linear solvers:
useBasis = or(opt.computeBasis, ~isempty(opt.mbasis));
[ls_pres, ls_TOF, ls_basis] = setupLinearSolvers(opt.linsolve, opt.reportTimings, useBasis);

% compute default basis if required:
basis = opt.mbasis;
if opt.computeBasis
    basis = computeDefaultBasis([], G, state, s, W, fluid, pv, T, 'linsolve', ls_basis);
    if opt.reportTimings
        t_basis = linsolveWithTimings;
    end
end

% evaluate initial state/objective
[states(1), Ds(1), grad] = solveStationaryPressure(G, state, s, W, fluid, pv, T, 'objective', objective,...
                    'linsolve', ls_pres, 'linsolveTOF', ls_TOF, 'msbasis', basis);
                
% set up objective scaling
sc = opt.objectiveScaling;
if isempty(sc) 
    v0 = grad.objective.val;
    sc = max(sqrt(eps), abs(v0));
end
scaling.boxLims = boxLims;
scaling.obj     = sc;     
% Set initial controls
uInit = well2control(W, 'scaling', scaling, 'targets', opt.targets);
% Scale linear equality cons
linEqSc = scaleCons(opt.linEq, scaling.boxLims);
opt.linEq = linEqSc;


% set up objective evalueation
f = @(u)evalObjectiveDiagnostics(u, objective, state, s, G, fluid, pv, T, W, scaling, ...
                                                'targets', opt.targets, ...
                                                'linSolve', ls_pres,    ...
                                                'linsolveTOF', ls_TOF,  ...
                                                'msbasis',  basis,      ...
                                                'minimize', opt.minimize);

[v, u, hst] = unitBoxBFGS(uInit, f, opt);

% final evaluation
[v, ~, W_opt, states(2), Ds(2)] = f(u);

disp(['Final objective value: ' num2str(v*sc)]);
info.history = hst;
info.scaling = scaling;
info.W_opt   = W_opt;
if opt.reportTimings
    t_lin = linsolveWithTimings;
    info.t_lin = reshape(t_lin, [6, numel(t_lin)/6]).';
    if opt.computeBasis
        info.t_basis = t_basis;
    end
end
end

function [ls_pres, ls_TOF, ls_basis] = setupLinearSolvers(linsolve, reportTimings, useBasis)
ls_basis = [];
if ~reportTimings
    if useBasis
        ls_pres  = @mldivide;
        ls_basis = linsolve;
    else
        ls_pres = linsolve;
    end
    ls_TOF  = @mldivide;
else
    trash = linsolveWithTimings; %#ok empty previously stored timings
    if useBasis
        ls_pres  = @linsolveWithTimings;
        ls_basis = @(A,x)linsolveWithTimings(A,x,linsolve);
    else
        ls_pres = @(A,x)linsolveWithTimings(A,x,linsolve); 
    end
    ls_TOF   = @linsolveWithTimings;
end
end

function linEqSc = scaleCons(linEq, bx)
if ~isempty(linEq)
    linEqSc.A = linEq.A*diag(bx(:,2)-bx(:,1));
    linEqSc.b = linEq.b - linEq.A*bx(:,1);
else
    linEqSc = [];
end
end

