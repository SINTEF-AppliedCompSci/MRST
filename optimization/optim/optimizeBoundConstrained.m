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
%                            Default: false
%       'stepInit'         : Initial step (gradient scaling). If not provided 
%                            (or set to nan), the following scaling will be used:
%                            stepInit = maxInitialUpdate/max(|initial gradient|)
%       'maxInitialUpdate' : as described above. Default: 0.05
%       'history'          : Option for 'warm starting' based on previous
%                            optimization (must have been run with option 
%                            'outputHessian' = true)
%       'lb'               : Lower bound(s) on u. If ~= 0, problem will be scaled
%                            Default: 0
%       'ub'               : Upper bound(s) on u. If ~= 1, problem will be scaled
%                            Default: 1
%   Stopping criteria options
%       'gradTol'          : Absolute tollerance of inf-norm of projected gradient. 
%                            Default: 1e-3
%       'objChangeTol'     : Absolute objective update tollerance. Default: 5e-4
%       'maxIt'            : Maximal number of iterations. Default: 25
%   Line-search options
%       'lineSearchMaxIt'  : Maximal number of line-search iterations. Default: 5
%       'stepIncreaseTol'  : Maximal step increase factor between line
%                            search iterations. Deafult: 10
%       'wolfe1',          : Objective improvement condition. Default: 1e-4
%       'wolfe2',          : Gradient reduction condition. Default: 0.9
%   Hessian approximation options
%       'lbfgsNum'         : Number of vector-pairs stored for L-BFGS.
%                          : Default: 10
%       'lbfgsStrategy'    : 'static' or 'dynamic' (see LimitedMemoryHessian.m)
%                            Default: 'dynamic'
%       'lbfgsRequireWolfe': Require that Wolfe conditions are satisfied
%                            in order to do BFGS-update. Default: false 
%   QP-solve options
%       'maxItQP',         : Maximal number of iterations for solving QP-problem
%                            for obtaining search direction. Related to how many 
%                            new bounds become active during an iteration.
%                            Deafult: 250
%       'activeChunkTol'   : Option to adjust for a more rough search reducing 
%                            the number its above. Default: sqrt(eps)         
%   Trust region options
%       'useTrustregion'   : Use infinity-norm trust region. Default: false
%       'radiusIncrease'   : Update-factor for 'good' approximations. Deafault: 2 
%       'radiusDecrease'   : Update-factor for 'bad' approximations. Default: 1/4
%       'ratioThresholds'  : rho-value thresholds for when radius will be updated
%                            i.e., decrease for < fac1, increase for > fac2. 
%                            Default: [.25 .75] 
%       'trustRegionInit'  : Initial trust-region radius. Default: 'stepInit'
%   Plotting and output options
%       'plotEvolution'    : Plot progess of optimization in figure. 
%                            Default: true
%       'logPlot'          : Use logarithmic y-axis when plotting objective
%                            Default: false
%       'outputHessian'    : Output Hessian approximation for each iteration
%                            in history-structure. Deafault: false 
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
%            'rho'   : local quadratic model fit
%            'r'     : trust region radius
%            'qpit'  : number of qp iterations
%            'qpit2' : number of outer qp iterations

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
                'lb',                  0,    ...
                'ub',                  1,    ...
                'stepInit',            nan,   ...
                'maxInitialUpdate',    0.05, ...
                'gradTol',             1e-3, ...
                'objChangeTol',        5e-4, ...
                'objChangeTolRel',     -inf, ...
                'maxIt',               25,   ...
                'lineSearchMaxIt',     5,   ...       
                'wolfe1',              1e-4, ...
                'wolfe2',              0.9,  ...
                'safeguardFac',        1e-5, ...
                'stepIncreaseTol',     10,   ...
                'maxItQP',             250,    ...
                'activeChunkTol',      sqrt(eps),   ... 
                'lbfgsNum',            10, ...
                'lbfgsStrategy',       'dynamic', ...
                'lbfgsRequireWolfe',   false, ...
                'useTrustRegion',      false, ... 
                'trustRegionInit',     nan, ...
                'radiusIncrease',      2,    ... 
                'radiusDecrease',      1/4,  ... 
                'ratioThresholds',     [.25 .75], ... 
                'plotEvolution',       true, ...
                'logPlot',             false, ... 
                'outputHessian',       false, ...
                'history',             [] );
opt  = merge_options(opt, varargin{:});
% negate f if we are maximizing
objSign = 1;
if opt.maximize
    f = @(u)fNegative(u, f);
    objSign = -1;
end
% scale problem to 0 <= u <= 1
objScale = false;
if any(opt.lb~= 0) || any(opt.ub~=1)
    f = @(u)fScale(u, f, opt.lb, opt.ub);
    objScale = true;
    u0 = (u0-opt.lb)./(opt.ub-opt.lb);
end
[lb, ub] = deal(0,1);

if isempty(opt.history) % starting from scratch
    if any(lb > u0) || any(u0 > ub)
        warning('Initial guess was not within bounds, projecting to feasible domain.')
        u0 = max(lb, min(ub, u0));
    end
    % Perform initial evaluation of objective and gradient:
    [v0, g0] = f(u0);
    % If not provided, set initial step
    step = opt.stepInit;
    if isnan(step) || step <= 0
        step = opt.maxInitialUpdate/max(abs(g0));
    end
    rTrust = opt.trustRegionInit;
    if opt.useTrustRegion && isnan(rTrust)
        rTrust = opt.maxInitialUpdate;
    end
    Hi = LimitedMemoryHessian('initScale', step, ...
                              'm', opt.lbfgsNum, ...
                              'initStrategy', opt.lbfgsStrategy);  
    HiPrev = Hi;
    it = 0;
    % Setup struct for gathering optimization history
    history = gatherInfo([], objSign*v0, u0, norm(g0), nan, nan, nan, Hi, nan, rTrust, nan, nan, opt.outputHessian);
else % starting from previous optimization
    history = opt.history;
    it = numel(history.val);
    Hi  = history.hess{it};
    assert(~isempty(Hi), ...
        'Warm start based on history requires Hessian approximations');
    u0 = history.u{it};
    [v0, g0] = f(u0);
    HiPrev = history.hess{max(it-1,1)};
    opt.maxIt = opt.maxIt + it;
    rTrust = history.r(it);
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
    if ~opt.useTrustRegion
        [lbcur, ubcur] = deal(lb, ub);
    else
        [lbcur, ubcur] = incorporateTrustRegion(u0, rTrust, lb, ub);
    end
    [d, Hi, pg, opt.maxStep, dObjEst, qpinfo] = getSearchDirection(u0, g0, Hi, HiPrev, lbcur, ubcur, opt);
    if ~(norm(pg,inf) < opt.gradTol) && ~isempty(d)
        [u, v, g, lsinfo] = lineSearch(u0, v0, g0, d, f, opt);
        dObjTrue = lsinfo.objVals(1) - v0;
        rho      = dObjTrue/dObjEst;
        if opt.useTrustRegion  
            rTrust = updateTrustRegion(rTrust, rho, norm(d, Inf), lsinfo.step, opt);
        end
        % check requirements for updating Hessian
        [du, dg] = deal(u-u0, g-g0);
        doUpdate = du'*dg > sqrt(eps)*norm(du)*norm(dg);
        if opt.lbfgsRequireWolfe
            doUpdate = doUpdate && lsinfo.flag > 0;
        end
        if  doUpdate
            % if any of the gradient entries are not defined, we set
            % the difference to zero
            dg(~isfinite(dg)) = 0;
            HiPrev = Hi;
            Hi     = Hi.update(du, dg);
        else
            fprintf('Hessian not updated during iteration %d.\n', it)
        end
        % Update history
        history = gatherInfo(history, objSign*v, u, norm(pg,inf), lsinfo.step, ...
            lsinfo.nits, lsinfo.flag, Hi, rho, rTrust, qpinfo.nit, qpinfo.nitOuter, opt.outputHessian);
    else
        if it == 1, [u,v, rho] = deal(u0, v0, nan); end
        history = gatherInfo(history, objSign*v, u, norm(pg,inf), 0, 0, 0, ...
            Hi, rho, rTrust, qpinfo.nit, qpinfo.nitOuter, opt.outputHessian);
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
end
if objScale
    u = u.*(opt.ub-opt.lb)+opt.lb;
end
end

%--------------------------------------------------------------------------
function [d, Hi, pg, maxStep, dObj, qpinfo] = getSearchDirection(u, g, Hi, HiPrev, lb, ub, opt)
% Find search-direction d which is the minimum of quadratic model 
%   m(d) = d'*g + d'*H*d, s.t lb <= u+d <= ub
for nTrial = 1:3
    if nTrial==2
        Hi = HiPrev;
    elseif nTrial==3
        Hi = Hi.reset();
    end
    % check whether projected gradient is below thresold
    active = getActiveBounds(u, -g, lb, ub);
    pg     = g(~active);
    if norm(pg, Inf) <= opt.gradTol
        [d, maxStep, dObj] = deal([], 0, 0);
        qpinfo = struct('nit', 0, 'nitOuter', 0);
        return;
    end
    
    [d, nit, nitOuter, dObj] = deal(0);
    done = false;
    gr   = g;    
    while ~done && nit < opt.maxItQP
        nitOuter = nitOuter + 1;
        active = getActiveBounds(u+d, -gr, lb, ub);
        % at minumum of the quadratic model, the gradient should be zero
        % except at active bounds
        [done, doneInner] = deal( norm(gr(~active), inf) < sqrt(eps) );
        while ~doneInner && nit < opt.maxItQP
            nit = nit +1;
            % get next (peace of) search direction from u+d
            dr  = -projQ(gr, active, Hi);
            % find next bound (or chunk of bounds) that become active
            [ix, s] = findNextBounds(u+d, dr, active, lb, ub, opt.activeChunkTol);
            % if s >= 1, there are no more bounds to violate
            doneInner = s >= 1;
            % for s > 1, step is dr, otherwise s*dr
            s   = min(1, s);
            sdr = s*dr;
            if numel(ix) > 1 && s > 0 && opt.activeChunkTol > 0
                % we may have slight (~tol) bound violations 
                sdr  = max(lb, min(ub, u+d+sdr)) - (u+sdr);
            end
            % update search direction:
            d = d + sdr;   
            % add quad model objective improvement for current segment:
            dObj = dObj + (1-s/2)*(sdr'*gr);
            % obtain gradient at u+d for quadratic model
            gr = gr - s*projQ(gr, active);
            % add new active bound(s) for next inner it if we are not done
            if ~doneInner
                active(ix) = true;
            end
        end
    end
    if ~done
        warning('Unable to solve local QP-problem in %d iterations.', opt.maxItQP);
    end
    d  = max(lb, min(ub, u+d)) - u;
    % find max step size before hitting next bound
    [~, maxStep] = findNextBounds(u, d, false(size(u)), lb, ub, 0);
    if maxStep < 1-sqrt(eps)
        warning('Problematic search direction, maximum step: %f < 1\n', maxStep);
    end
    qpinfo = struct('nit', nit, 'nitOuter', nitOuter);
    isDecreasing  = d'*g <= 0;
    if  isDecreasing % decreasing search direction found, exit
        break;
    else
        % retry with other Hessian approx
        str = 'Non-inceasing search direction';
        switch nTrial
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

function w = projQ(v, active, Hi)
% computer action of inverse of Hessian restricted to non-active
if nargin < 3
    w = v.*(~active);
else
    if any(active)
        Hi = Hi.setNullspace(active);
    else
        Hi = Hi.setNullspace();
    end
    w = Hi*v;
end
end
    
%--------------------------------------------------------------------------
function [active, type] = getActiveBounds(u, v, lb, ub)
% get active bound-flag for each variable
activeLower = (u < lb + sqrt(eps)) & (v < 0);
activeUpper = (u > ub - sqrt(eps)) & (v > 0);
active = activeUpper | activeLower;
if nargout == 2
    type = -activeLower + activeUpper;
end
end

%--------------------------------------------------------------------------
function [ix, smax] = findNextBounds(u, d, active, lb, ub, tol)
% Find next occuring bound(s) from u in direction d. d is zero for
% active variables, but might also be zero for others
[sl, su] = deal((lb-u)./d, (ub-u)./d);
s  = max(sl, su); % pick whichever is positive
% disregard d == 0
s(active | d == 0) = inf;
[smax, ix] = min(s);
if tol > 0 && smax < 1
   % find all within tollerance
   ix = find(s <= smax*(1+tol));
   % select maximum (i.e., all become active/violated)
   smax = max(s(ix));
end
end


%--------------------------------------------------------------------------
function [v, g] = fNegative(u, f)
% negate function
if nargout == 1
    v = -f(u);
elseif nargout == 2
    [v, g] = f(u);
    [v, g] = deal(-v, -g);
end
end

%--------------------------------------------------------------------------
function [v, g] = fScale(u, f, lb, ub)
% scale function
if nargout == 1
    v = f( u.*(ub-lb)+lb );
elseif nargout == 2
    [v, g] = f( u.*(ub-lb)+lb );
    g = g.*(ub-lb);
end
end

%--------------------------------------------------------------------------
function [lbcur, ubcur] = incorporateTrustRegion(u, rTrust, lb, ub)
% update bounds for current trust region radius (inf-norm)
lbcur = max(lb, u-rTrust);
ubcur = min(ub, u+rTrust);
end

%--------------------------------------------------------------------------
function r = updateTrustRegion(r, rho, update, step, opt)
% update trust-region based on rho (quad-model fit)
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
function hst = gatherInfo(hst, val, u, pg, alpha, lsit, lsfl, hess, rho, r, qpit, qpit2, outputH)   
% obj.val | contr | norm proj.grad | ls-step | ls-its | ls-flag | hessian 
if ~outputH
    hess = [];
end
if isempty(hst)
    hst = struct('val', val, 'u', {{u}}, 'pg', pg, 'alpha', alpha, ...
                 'lsit', lsit, 'lsfl', lsfl, 'hess', {{hess}}, ...
                 'rho', rho, 'r', r, 'qpit', qpit, 'qpit2', qpit2);
else
    hst.val   = [hst.val  , val  ];
    hst.u     = [hst.u    , {u}  ];
    hst.pg    = [hst.pg   , pg   ];
    hst.alpha = [hst.alpha, alpha];
    hst.lsit  = [hst.lsit , lsit ];
    hst.lsfl  = [hst.lsfl , lsfl ];
    hst.hess  = [hst.hess ,{hess}];
    hst.rho   = [hst.rho  , rho  ];
    hst.r     = [hst.r, r];
    hst.qpit  = [hst.qpit, qpit];
    hst.qpit2 = [hst.qpit2, qpit2];
end
end

%--------------------------------------------------------------------------
function printInfo(history, it)
fprintf('It: %2.1d | val: %4.3e | ls-its: %3.1d | pgrad: %4.3e | rho: %4.3e\n', ...
        it, history.val(end), history.lsit(end), ...
        history.pg(end), history.rho(end));
end
