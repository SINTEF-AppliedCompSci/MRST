function [v, u, history] = unitBoxBFGS(u0, f, varargin)
% Iterative line search optimization using BFGS intended for scaled 
% problems, i.e., 0<=u<=1 and f~O(1) 
% 
% SYNOPSIS:
% [v, u, history] = unitBoxBFGS(u0, f, ...)
%
% PARAMETERS
% u0    : inital guess nx1 vector with 0<= u0 <= 1 and feasible wrt
%         any additional constraints 
% f     : function handle s.t., [v,g] = f(u) returns
%           v : objective value
%           g : objective nx1 gradient vector
% KEYWORD ARGUMENTS:
%       'maximize'         : Boolean option (false will minimize objective)
%                            Default: true
%       'stepInit'         : Initial step (gradient scaling). If not provided 
%                            (or set to nan), the following scaling will be used:
%                            stepInit = maxInitialUpdate/max(|initial gradient|)
%       'maxInitialUpdate' : as described above. Default: 0.05
%       'history'          : Option for 'warm starting' based on previous
%                            optimization (must have been run with option 
%                            'outputHessian' = true)
%   Stopping criteria options
%       'gradTol'       : Absolute tollerance of inf-norm of projected gradient. 
%                         Default: 1e-3
%       'objChangeTol'  : Absolute objective update tollerance. Default: 5e-4
%       'maxIt'         : Maximal number of iterations. Default: 25
%   Line-search options
%       'lineSearchMaxIt' : Maximal number of line-search iterations. Default: 5
%       'stepIncreaseTol' : Maximal step increase factor between line
%                           search iterations. Deafult: 10
%       'wolfe1',         : Objective improvement condition. Default: 1e-4
%       'wolfe2',         : Gradient reduction condition. Default: 0.9
%   Hessian approximation options
%       'useBFGS'        : If false pure gradient search will be used.
%                          Deafult: true
%       'limitedMemory'  : If false, full Hessian approximations will be
%                          made (in general not recommended), otherwise
%                          L-BFGS (see LimitedMemoryHessian.m).
%                          Default: true
%       'lbfgsNum'       : Number of vector-pairs stored for L-BFGS.
%                        : Default: 5
%       'lbfgsStrategy'  : 'static' or 'dynamic' (see LimitedMemoryHessian.m)
%                          Default: 'dynamic'
%   Linear constraints options
%       'linEq'           : Linear equality constraints given as structure with
%                           fields 'A' and 'b' to enforce A*u=b
%       'linIneq'         : Linear inequality constraints *in addition to* the default 
%                           box constraints (0 <= u <= 1). Given as structure with
%                           fields 'A' and 'b' to enforce A*u<=b
%       'enforceFeasible' : If true, will attempt to repair problematic 
%                           constraint handling/violations by projection. A
%                           NOTE: A non-feasible initial guess, will always be
%                           attempted to be made feasible.
%                           Default: false
%   Plotting and output options
%       'plotEvolution' : Plot progess of optimization in figure. 
%                         Default: true
%       'logPlot'       : Use logarithmic y-axis when plotting objective
%                         Default: false
%       'outputHessian' : Output Hessian approximation for each iteration
%                         in history-structure. Deafault: false 
% 
% RETURNS:                          
% v       : Optimal or best objective value
% u       : Control/parameter vector corresponding to v
% history : Structure containing for each iteration:
%            'val'   : objective value
%            'u'     : control/parameter vector
%            'pg'    : norm of projected gradient
%            'alpha' : line-search step length
%            'lsit'  : number of line-search iterations
%            'lsfl'  : line-search flag
%            'hess'  : Hessian inverse approximation (if requested)

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
opt = struct(   'maximize',            true, ...
                'stepInit',            nan,   ...
                'maxInitialUpdate',    0.05, ...
                'gradTol',             1e-3, ...
                'objChangeTol',        5e-4, ...
                'maxIt',               25,   ...
                'lineSearchMaxIt',     5,   ...       
                'wolfe1',              1e-4, ...
                'wolfe2',              0.9,  ...
                'safeguardFac',        1e-5, ...
                'stepIncreaseTol',     10,   ...
                'useBFGS',             true, ...
                'limitedMemory',       true, ...
                'lbfgsNum',            5, ...
                'lbfgsStrategy',       'dynamic', ...
                'linEq',               [], ...
                'enforceFeasible',     true, ... 
                'linIneq',             [], ...
                'plotEvolution',       true, ...
                'logPlot',             false, ... 
                'outputHessian',       false, ...
                'history',             [] );
opt  = merge_options(opt, varargin{:});
objSign = 1;
if opt.maximize
    f = @(u)fNegative(u, f);
    objSign = -1;
end
c  = getConstraints(u0, opt);
if isempty(opt.history)
    % Setup constraint struct
    [u0, ~, consOK] = checkFeasible(u0, c, opt.enforceFeasible, 'Initial guess');
    assert(consOK, 'Infeasible initial guess')
    
    % Perform initial evaluation of objective and gradient:
    [v0,g0] = f(u0);
    [v ,u ] = deal(v0,u0);
    % If not provided, set initial step
    step = opt.stepInit;
    if isnan(step) || step <= 0
        step = opt.maxInitialUpdate/max(abs(g0));
    end
    % Initial Hessian-approximation
    if ~opt.limitedMemory
        Hi = step;
    else
        Hi = LimitedMemoryHessian('initScale', step, ...
            'm', opt.lbfgsNum, ...
            'initStrategy', opt.lbfgsStrategy);
    end
    HiPrev = Hi;
    it = 0;
    % Setup struct for gathering optimization history
    history = gatherInfo([], objSign*v0, u0, norm(g0), nan, nan, nan, Hi, opt.outputHessian);
else
    history = opt.history;
    it = numel(history.val);
    Hi  = history.hess{it};
    assert(~isempty(Hi), ...
        'Warm start based on history requires Hessian approximations');
    u0 = history.u{it};
    [v0, g0] = f(u0);
    HiPrev = history.hess{max(it-1,1)};
    opt.maxIt = opt.maxIt + it;
end
if opt.plotEvolution >0
    evFig = figure;
    plotInfo(evFig, history, opt.logPlot)
end
printInfo(history, it);

success = false;
g = g0;
while ~success
    it = it+1;
    [d, Hi, pg, opt.maxStep] = getSearchDirection(u0, g0, Hi, HiPrev, c);
    %if isempty(d), break; end
    if ~(norm(pg,inf) < opt.gradTol) && ~isempty(d)
        % Perform line-search
        [tmp, flag, fixed]  = checkFeasible(u0+d, c, opt.enforceFeasible);
        if ~flag && fixed
            d = tmp - u0;
        end
        [u, v, g, lsinfo] = lineSearch(u0, v0, g0, d, f, opt);
        % Update Hessian approximation
        if opt.useBFGS %&& lsinfo.flag == 1
            [du, dg] = deal(u-u0, g-g0);
            if du'*dg > sqrt(eps)*norm(du)*norm(dg)
                HiPrev = Hi;
                if isa(Hi, 'LimitedMemoryHessian')
                    Hi = Hi.update(du, dg);
                else
                    r = 1/(du'*dg);
                    V = eye(numel(u)) - r*dg*du';
                    Hi = V'*Hi*V + r*(du*du');
                end
            else
                fprintf('Hessian not updated during iteration %d.\n', it)
            end
        end
        % Update history
        history = gatherInfo(history, objSign*v, u, norm(pg,inf), lsinfo.step, ...
                             lsinfo.nits, lsinfo.flag, Hi, opt.outputHessian);
    else
        history = gatherInfo(history, objSign*v, u, norm(pg,inf), 0, 0, 0, ...
                             Hi, opt.outputHessian);
    end
    
    %Check stopping criteria
    success = (it >= opt.maxIt) || (norm(pg,inf) < opt.gradTol) || ...
              (abs(v-v0) < opt.objChangeTol);

    [u0, v0, g0] = deal(u, v, g);
    if opt.plotEvolution > 0
        plotInfo(evFig, history, opt.logPlot)
    end
    printInfo(history, it);
end
end
%--------------------------------------------------------------------------

function [d, Hi, pg, maxStep] = getSearchDirection(u0, g0, Hi, HiPrev, c)
% Find search-direction which is (sum of) the projection(s) of Hi*g0 
% restricted to directions with non-active constraints. 
cnt = 1;
for k = 1:3
    if k==2
        Hi = HiPrev;
    elseif k==3
        if isa(Hi, 'LimitedMemoryHessian')
            Hi = Hi.reset();
        else
            Hi = 1;
        end
    end
    % Project gradient and search direction onto nullspace of equality 
    % constraints
    Q = c.e.Q;
    pg = -projQ(g0, Q);
    d  = -projQ(g0, Q, Hi);
    % The following loop(s) should project onto nullspace of additional 
    % active constraints. First for gradient direction, then for search 
    % direction which may have additional constraints active. While-loop 
    % will in most cases only be traversed once, so we don't worry too much 
    % about efficiency
    isActive = false;
    for kd = 1:2
        [na, na_prev] = deal(0, -1);
        while na > na_prev
            if kd == 1
                [~, active_cur] = classifyConstraints(c.i.A, c.i.b, u0, pg);
            else
                [~, active_cur] = classifyConstraints(c.i.A, c.i.b, u0, d);
            end
            isActive = isActive | active_cur;
            na_prev  = na;
            na       = nnz(isActive);
            if na > na_prev % redo projection for all active 
                [Q, s] = svd([c.i.A(isActive,:)', c.e.A'], 0);
                s  = diag(s);
                Q  = Q(:, s > sqrt(eps)*s(1));
                if kd == 1 
                    pg = -projQ(g0, Q);
                else       
                    d  = -projQ(g0, Q, Hi);
                end
            end
        end
    end
    % Check for tiny projected gradient
    if norm(pg, inf) <= sqrt(eps)*norm(g0, inf) %  nothing more to do
        [d, maxStep] = deal([]);
        return
    end
    % Iteratively find all constraints that become active from u0 to u0+d,
    % and project remaining line segments accordingly.
    [dr, gr] = deal(d, g0);
    becomesActive = isActive;
    d    = 0;
    done = false;
    while ~done
        if norm(dr) > sqrt(eps)
            sgn     = classifyConstraints(c.i.A, c.i.b, u0+d, dr);
            [ix, s] = findNextConstraint(c.i.A, c.i.b, u0+d, dr, sgn<=0 | becomesActive);
        else
            [ix, s] = deal([], 0);
        end
        if ~isempty(ix) && s <= 1+sqrt(eps)
            becomesActive(ix) = true;
            d  = d + s*dr;
            gr = (1-s)*gr;
            Q  = expandQ(Q, c.i.A(ix,:)');
            dr = -projQ(gr, Q, Hi);
        else
            d    = d+dr;
            done = true;
        end
    end
    % find maximal step length we can take with d before hitting the next
    % constraint (should be >= 1)
    sgn           = classifyConstraints(c.i.A, c.i.b, u0, d);
    [~, maxStep]  = findNextConstraint(c.i.A, c.i.b, u0, d, sgn<=0);
    if maxStep < .95
        warning('Problematic constraint handling, relative step length: %6.5f', maxStep);
    end
    if maxStep < 1
        [d, maxStep] = deal(maxStep*d, 1);
    end
   
    isDecreasing  = d'*g0 <= 0;
    isZero        = norm(d, inf) <= sqrt(eps)*norm(Hi*g0, inf);
    if  isDecreasing && ~isZero % decreasing search direction found, exit
        break;
    else
        % retry with other Hessian approx
        if ~isZero
            str = 'Small norm of search direction';
        else
            str = 'Non-inceasing search direction';
        end
        switch cnt
            case 1
                fprintf('%s, trying previous Hessian approximation.\n', str)
            case 2
                fprintf('%s, trying to reset Hessian to identity.\n', str)
            case 3
                fprintf('Exiting: %s.\n', str)
                [d, maxStep] = deal([]);
        end
    end
end
end
%--------------------------------------------------------------------------

function w = projQ(v, Q, H)
if nargin < 3
    H = 1;
end
if isempty(Q)
    w = H*v;
else
    if ~isa(H, 'LimitedMemoryHessian')
        tmp = H*(v-Q*(Q'*v));
        w   = tmp - Q*(Q'*tmp);
    else
        H = H.setNullspace(Q);
        w = H*v;
    end
end
end
%--------------------------------------------------------------------------

function c = getConstraints(u, opt)
% Box constraints, always 0<= u <= 1
nu = numel(u);
c.i.A = [diag(-ones(nu,1)); diag(ones(nu,1))];
c.i.b = [zeros(nu,1); ones(nu,1)];
% Add general linear constraints 
if ~isempty(opt.linIneq)
    % scale by norm(A)
    sc = norm(opt.linIneq.A);
    c.i.A = [c.i.A; opt.linIneq.A/sc];
    c.i.b = [c.i.b; opt.linIneq.b/sc];
end
% Equality constraints (always active)
if ~isempty(opt.linEq)
    sc = norm(opt.linEq.A);
    c.e.A = opt.linEq.A/sc;
    c.e.b = opt.linEq.b/sc;
    [Q, s] = svd(c.e.A', 0);
    s      = diag(s);
    c.e.Q  = Q(:, s> sqrt(eps) * s(1));
else
    c.e.A = zeros(0, nu);
    c.e.b = [];
    c.e.Q = zeros(nu, 0);
end
end
%--------------------------------------------------------------------------

function [u, flag, fixed] = checkFeasible(u, c, enforce, nm)
% Check that u is feasible. If not and enforce == true, try to fix. Ideally
% should be solved as ||u-u*|| st c.e and c.i, but since we don't have a
% QP-solver resort to (costly) iterative projections in active subspaces. 
% Intended for fixing mild violations.
if nargin < 4
    nm = 'Vector u';
end
if nargin < 3
    enforce = false;
end
hasEC = ~isempty(c.e.A);
hasIC = ~isempty(c.i.A);
[ecOK, icOK] = deal(true);

if hasEC
    [Ae, be] = deal(c.e.A, c.e.b);
    if any(abs(c.e.A*u-c.e.b) > sqrt(eps))
        % find closest u that fulfills c.e.
        u = u + Ae'*((Ae*Ae')\(be-Ae*u));
        ecOK = false; %now OK, but warn
    end
end

flag  = ecOK;
fixed = false;
maxIt = 100;
for it = 1:maxIt
    if hasIC
        icOK = ~any(c.i.A*u-c.i.b > sqrt(eps));
        flag = flag && icOK;
    end
    if ~enforce
        break
    else
        if ~icOK
            Q    = zeros(numel(u), 0);
            proj = @(v)v;
            if hasEC
                Q     = c.e.Q;
                proj  = @(v)v-Q*(Q'*v);
            end
            % Loop through each violating and project. If there are no more
            % availabe directions, restart while-loop from current point
            done = false;
            cnt  = 0;
            while ~done
                if cnt == 0
                    icIx = find(c.i.A*u-c.i.b > sqrt(eps));
                    icIx = circshift(icIx, it);
                end
                if isempty(icIx) || cnt == numel(icIx)
                    done = true;
                else
                    %
                    ix   = icIx(cnt+1);
                    [a, b]  = deal(c.i.A(ix,:)', c.i.b(ix));
                    pa = proj(a);
                    if norm(pa) < sqrt(eps)*norm(a)
                        % skip
                        cnt = cnt +1;
                    else
                        cnt = 0;
                        u = u + pa*((b-a'*u)/(a'*pa));
                        if size(Q,2) < size(Q,1)-1
                            Q = expandQ(Q, pa);
                            proj = @(v)v-Q*(Q'*v);
                        else
                            % possibly restart while-loop
                            done = true;
                        end
                    end
                end
            end
        else
            fixed = true;
            break
        end
    end
end

 
if it == maxIt
    warning('Failed attempt to fix feasibility of %s, continuing anyway ...', nm);
elseif ~flag
    if ~enforce
        warning('%s is not feasible within tollerance. %s', nm, ...
                'Consider running with option ''enforceFeasible''=true');
    else
        warning('%s was not feasible, fixed feasibility in %d iteration(s)', ...
                nm, it-1);
    end
end
end
%--------------------------------------------------------------------------
function [sgn, act] = classifyConstraints(A, b, u, v)
% classify inequality constraints for point u with direction v
% sgn: -1: in, 0: parallell, 1: out
% act: true for active
sgn = A*v;
sgn(abs(sgn)<sqrt(eps)) = 0;
sgn = sign(sgn);
act = A*u-b > -sqrt(eps) & sgn > 0;
end
%--------------------------------------------------------------------------

function [ix, s] = findNextConstraint(A, b, u, d, ac)
s = (b-A*u)./(A*d);
s(ac)  = inf;
s(s<eps) = inf;
[s, ix] = min(s);
end
%--------------------------------------------------------------------------

function Q = expandQ(Q, v)
n0 = norm(v);
v = v - Q*(Q'*v);
if norm(v)/n0 > sqrt(eps)
    Q = [Q, v/norm(v)];
else
    fprintf('Newly active constraint is linear combination of other active constraints ??!!\n')
end
end
%--------------------------------------------------------------------------

function [v, g] = fNegative(u, f)
if nargout == 1
    v = -f(u);
elseif nargout == 2
    [v, g] = f(u);
    [v, g] = deal(-v, -g);
end
end
%--------------------------------------------------------------------------

function [] = plotInfo(fig, hst, logPlot)
if ~ishandle(fig)
    figure(fig)
else
    % Avoid stealing focus if figure already exists
    set(0, 'CurrentFigure', fig);
end
xt = 0:(numel(hst.val)-1);
xlim = [-.2, xt(end)+.5];
ch = abs(hst.val(2:end)-hst.val(1:end-1));
popt = {'o-', 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerFaceColor', [1 1 1]};
subplot(5,1,1), 
if ~logPlot
    plot(xt, hst.val, popt{:});
else
    semilogy(xt, abs(hst.val), popt{:});
end
title('Objective');
set(gca, 'XLim', xlim)
subplot(5,1,2), semilogy(xt,hst.pg, popt{:}), title('Gradient norm');
set(gca, 'XLim', xlim)
subplot(5,1,3), semilogy(xt(2:end)-.5,ch, popt{:}), title('Objective change');
set(gca, 'XLim', xlim)
subplot(5,1,4), bar(xt,hst.lsit), title('Line search iterations');
set(gca, 'XLim', xlim)
subplot(5,1,5), col = [0.8500 0.3250 0.0980];
if numel(hst.u{end}) < 50
    bar(hst.u{end}, 'FaceColor', col);
else
    plot(hst.u{end}, '.-', 'Color', col);
end
title('Current scaled controls');
set(gca, 'YLim', [0, 1])
drawnow
end
%--------------------------------------------------------------------------

function hst = gatherInfo(hst, val, u, pg, alpha, lsit, lsfl, hess, outputH)   
% obj.val | contr | norm proj.grad | ls-step | ls-its | ls-flag | hessian 
if ~outputH
    hess = [];
end
if isempty(hst)
    hst = struct('val', val, 'u', {{u}}, 'pg', pg, ...
                 'alpha', alpha, 'lsit', lsit, 'lsfl', lsfl, ...
                 'hess', {{hess}});
else
    hst.val   = [hst.val  , val  ];
    hst.u     = [hst.u    , {u}  ];
    hst.pg    = [hst.pg   , pg   ];
    hst.alpha = [hst.alpha, alpha];
    hst.lsit  = [hst.lsit , lsit ];
    hst.lsfl  = [hst.lsfl , lsfl ];
    hst.hess  = [hst.hess ,{hess}];
end
end
%--------------------------------------------------------------------------

function printInfo(history, it)
fprintf('It: %2.1d | val: %4.3e | ls-its: %3.1d | pgrad: %4.3e\n', ...
        it, history.val(end), history.lsit(end), ...
        history.pg(end));
end
