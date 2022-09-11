function [v, u, history] = optimizeBoundConstrained(u0, f, varargin)
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
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
                'trustRegion',         true, ... 
                'lb',                  0,    ...
                'ub',                  1,    ...
                'stepInit',            nan,   ...
                'maxInitialUpdate',    0.05, ...
                'gradTol',             1e-8, ...
                'objChangeTol',        5e-8, ...
                'objChangeTolRel',     -inf, ...
                'maxIt',               25,   ...
                'lineSearchMaxIt',     5,   ...       
                'wolfe1',              1e-4, ...
                'wolfe2',              0.9,  ...
                'safeguardFac',        1e-5, ...
                'stepIncreaseTol',     10,   ...
                'radiusIncrease',      2,    ... % increase radius  (good approx)
                'radiusDecrease',      1/4,  ... % decrease radius  (bad approx)
                'ratioThresholds',     [.25 .75], ... % bad | medium | good approx
                'maxItQP',             100,    ...
                'activeChunkTol',      1e-3,   ... 
                'lbfgsNum',            5, ...
                'lbfgsStrategy',       'dynamic', ...
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
[lb, ub] = deal(opt.lb, opt.ub);
objScale = false;
if any(lb~= 0) || any(ub~=1)
    f = @(u)fScale(u, f, lb, ub);
    objScale = true;
end
if isempty(opt.history)
    if any(lb > u0) || any(u0 > ub)
        warning('Initial guess was not within bounds, projecting to feasible domain.')
        u0 = max(lb, min(ub, u0));
    end
    % Perform initial evaluation of objective and gradient:
    [v0,g0] = f(u0);
    % If not provided, set initial step
    step = opt.stepInit;
    if isnan(step) || step <= 0
        step = opt.maxInitialUpdate/max(abs(g0));
    end
    rTrust = step;
    Hi = LimitedMemoryHessian('initScale', step, ...
                              'm', opt.lbfgsNum, ...
                              'initStrategy', opt.lbfgsStrategy);
  
    HiPrev = Hi;
    it = 0;
    % Setup struct for gathering optimization history
    history = gatherInfo([], objSign*v0, u0, norm(g0), nan, nan, nan, Hi, nan, opt.outputHessian);
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
[v, u] = deal(v0,u0);
if opt.plotEvolution >0
    evFig = figure;
    plotInfo(evFig, history, opt.logPlot)
end
printInfo(history, it);

success = false;
g = g0;
while ~success
    it = it+1;
    if ~opt.trustRegion
        [lbcur, ubcur] = deal(lb, ub);
    else
        [lbcur, ubcur] = incorporateTrustRegion(u0, rTrust, lb, ub);
    end
    [d, Hi, pg, opt.maxStep, dObjEst] = getSearchDirection(u0, g0, Hi, HiPrev, lbcur, ubcur, opt);
    if ~(norm(pg,inf) < opt.gradTol) && ~isempty(d)
        [u, v, g, lsinfo] = lineSearch(u0, v0, g0, d, f, opt);
        dObjTrue = lsinfo.objVals(1) - v0;
        rho      = dObjTrue/dObjEst;
        if opt.trustRegion  
            rTrust = updateTrustRegion(rTrust, rho, norm(d, Inf), lsinfo.step, opt);
        end
        % Update Hessian approximation
        [du, dg] = deal(u-u0, g-g0);
        % if any of the gradient entries are not defined, we set
        % the difference to zero
        dg(~isfinite(dg)) = 0;
        if du'*dg > sqrt(eps)*norm(du)*norm(dg) %&& lsinfo.flag > 0
            HiPrev = Hi;
            Hi = Hi.update(du, dg);
        else
            fprintf('Hessian not updated during iteration %d.\n', it)
        end
        % Update history
        history = gatherInfo(history, objSign*v, u, norm(pg,inf), lsinfo.step, ...
            lsinfo.nits, lsinfo.flag, Hi, rho, opt.outputHessian);
        %        history = gatherInfo(history, objSign*v, u, norm(pg,inf), rho, ...
        %    1, 1, Hi, opt.outputHessian);
    else
        history = gatherInfo(history, objSign*v, u, norm(pg,inf), 0, 0, 0, ...
                             Hi, rho, opt.outputHessian);
    end
    
    %Check stopping criteria
    success = (it >= opt.maxIt) || (norm(pg,inf) < opt.gradTol) || ...
              (abs(v-v0) < opt.objChangeTol) || ...
              (abs((v-v0)./v) < opt.objChangeTolRel);

    [u0, v0, g0] = deal(u, v, g);
    if opt.plotEvolution > 0
        plotInfo(evFig, history, opt.logPlot)
    end
    printInfo(history, it);
    if objScale
        u = u*(ub-lb)+lb;
    end
end
end
%--------------------------------------------------------------------------

function [d, Hi, pg, maxStep, dObj] = getSearchDirection(u, g, Hi, HiPrev, lb, ub, opt)
% Find search-direction which is (sum of) the projection(s) of Hi*g0 
% restricted to directions with non-active constraints. 
cnt = 1;
for k = 1:3
    if k==2
        Hi = HiPrev;
    elseif k==3
            Hi = Hi.reset();
    end
    d = 0;
    % output projected gradient at current point
    active = getActiveBounds(u+d, -g, lb, ub);  
    pg = projQ(g, active);
    if norm(pg, Inf) <= opt.gradTol
        [d, maxStep, dObj] = deal(nan, 0, 0);
        return;
    end
    [done, nit] = deal(false, 0);
    dObj = 0; % predicted improvement in objective
    gr = g;
    active = false(size(u));
    d0 = -projQ(g, active, Hi); % main search direction
    while ~done && nit < opt.maxItQP
        nit = nit +1;
        activePrev = active;
        active = getActiveBounds(u+d, d0, lb, ub);
        dr     = -projQ(gr, active, Hi);
        [ix, smax] = findNextBounds(u+d, dr, active, lb, ub, opt.activeChunkTol);
        if  ~isempty(ix) && smax < 1
            %active(ix) = true;
            d  = d + smax*dr;
            % estimated improvement
            dObj = dObj + smax*(1-.5*smax)*(dr'*gr);
            % update remaining gradient/step
            gr = gr - smax*projQ(gr, active);
        else
            d = d + dr;
            dObj = dObj + .5*(dr'*gr);
            gr = gr.*active;
            if all(active == activePrev)
                done = true;
            end
        end
    end
    if ~done
        warning('Bla');
    end
    d  = max(lb, min(ub, u+d)) - u;
    % find max step size before hitting next bound
    [~, maxStep] = findNextBounds(u, d, false(size(u)), lb, ub, 0);
    if maxStep < 1-sqrt(eps)
        warning('Problematic search direction, debug ...');
    end
    isDecreasing  = d'*g <= 0;
    if  isDecreasing % decreasing search direction found, exit
        break;
    else
        % retry with other Hessian approx
        str = 'Non-inceasing search direction';
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

function w = projQ(v, active, H)
if nargin < 3
    w = v.*(~active);
else
    if any(active)
        H = H.setNullspace(active);
    else
        H = H.setNullspace();
    end
    w = H*v;
end
end
    
%--------------------------------------------------------------------------
function [active, type] = getActiveBounds(u, v, lb, ub)
activeLower = (u < lb + sqrt(eps)) & (v < 0);
activeUpper = (u > ub - sqrt(eps)) & (v > 0);
active = activeUpper | activeLower;
if nargout == 2
    type = -activeLower + activeUpper;
end
end

%--------------------------------------------------------------------------
function [ix, smax] = findNextBounds(u, d, active, lb, ub, tol)
d(d==0) = eps;
[sl, su] = deal((lb-u)./d, (ub-u)./d);
s  = max(sl, su); % pick whichever is positive
s(active) = inf;
[smax, ix] = min(s);
if tol > 0 && smax < 1
   ix = find(s <= smax*(1+tol));
   smax = max(s(ix));
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

function [v, g] = fScale(u, f, lb, ub)
if nargout == 1
    v = f( u.*(ub-lb)+lb );
elseif nargout == 2
    [v, g] = f( u.*(ub-lb)+lb );
    g = g.*(ub-lb);
end
end

%--------------------------------------------------------------------------
function [lbcur, ubcur] = incorporateTrustRegion(u, rTrust, lb, ub)
lbcur = max(lb, u-rTrust);
ubcur = min(ub, u+rTrust);
end
%--------------------------------------------------------------------------
function r = updateTrustRegion(r, rho, update, step, opt)
if rho < opt.ratioThresholds(1)
    r = opt.radiusDecrease*step*r;
elseif rho > opt.ratioThresholds(2) && r < update*step*(1+sqrt(eps))
    r = opt.radiusIncrease*step*r;
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

function hst = gatherInfo(hst, val, u, pg, alpha, lsit, lsfl, hess, rho, outputH)   
% obj.val | contr | norm proj.grad | ls-step | ls-its | ls-flag | hessian 
if ~outputH
    hess = [];
end
if isempty(hst)
    hst = struct('val', val, 'u', {{u}}, 'pg', pg, ...
                 'alpha', alpha, 'lsit', lsit, 'lsfl', lsfl, ...
                 'hess', {{hess}}, 'rho', rho);
else
    hst.val   = [hst.val  , val  ];
    hst.u     = [hst.u    , {u}  ];
    hst.pg    = [hst.pg   , pg   ];
    hst.alpha = [hst.alpha, alpha];
    hst.lsit  = [hst.lsit , lsit ];
    hst.lsfl  = [hst.lsfl , lsfl ];
    hst.hess  = [hst.hess ,{hess}];
    hst.rho   = [hst.rho  , rho  ];
end
end
%--------------------------------------------------------------------------

function printInfo(history, it)
fprintf('It: %2.1d | val: %4.3e | ls-its: %3.1d | pgrad: %4.3e | rho: %4.3e\n', ...
        it, history.val(end), history.lsit(end), ...
        history.pg(end), history.rho(end));
end
