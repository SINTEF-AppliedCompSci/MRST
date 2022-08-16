function [v, u, history] = unitBoxLM(u0, f, varargin)
% Levenberg-Marquardt-type optimization intended for scaled 
% problems, i.e., 0<=u<=1 and f~O(1) 
% SYNOPSIS:
% [v, u, history] = unitBoxLM(u0, f, ...)
%
% DESCRIPTION:
% See options. There are two options for lambda-updating:
% (1)'updateStrategy' = 'simple' (default) will perform a quasi-trust-region
% update directly on lambda using 'lambdaIncrease'/'lambdaDecrease'
% (2)'updateStrategy' = 'TR' will perform a trust-region update on norm of
% update using 'radiusIncrease'/'radiusDecrease' and compute an
% approximative corresponding lambda-value.
%
% 
% PARAMETERS
% u0    : inital guess nx1 vector with 0<= u0 <= 1 
% f     : function handle s.t., [v,J] = f(u) returns
%           v : nx1 vector of residuals
%           J : nxm Jacobian such that J'*v is gradient of v.^2 wrt u
% KEYWORD ARGUMENTS:
%        See below
% RETURNS:                          
% v       : Optimal or best objective value
% u       : Control/parameter vector corresponding to v
% history : Structure containing for each iteration:
%            'val'    : objective value (sum(v.^2))
%            'u'      : control/parameter vector
%            'pg'     : norm of projected gradient
%            'lambda' : damping parameter
%            'rho'    : local model agreement
%            'nIt'    : number of local iterations

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
opt = struct('lambdaInit',               0.01, ...
             'lambdaIncrease',              8, ... % increase damping (bad approx)
             'lambdaDecrease',              5, ... % decrease damping (good approx)
             'radiusIncrease',              2, ... % increase radius  (good approx)
             'radiusDecrease',              4, ... % decrease radius  (bad approx)
             'lambdaMax',                 1e6, ...
             'lambdaMin',                1e-6, ...
             'ratioThresholds',     [.25 .75], ... % bad | medium | good approx
             'scaledDamping',           false, ... % not recommended
             'updateStrategy',       'simple', ... % 'simple' or 'TR' (trust region)
             'gradTol',                  1e-6, ... 
             'updateTol',                1e-6, ...
             'resTolAbs',                1e-5, ...
             'resTolRel',                   0, ...
             'resChangeTolRel',          -inf, ...
             'maxIt',                      20, ...
             'maxFunEvals',                [], ...
             'lsqTol',                      0, ... % if >0, use lsqminnorm
             'plotEvolution',            true, ...
             'verbose',                  true);
 
opt  = merge_options(opt, varargin{:});
if isempty(opt.maxFunEvals)
    opt.maxFunEvals = 2*opt.maxIt;
end
if ~strcmp(opt.updateStrategy, 'TR')
    goodFac = opt.lambdaDecrease;
    badFac  = opt.lambdaIncrease;
else
    goodFac = opt.radiusIncrease;
    badFac  = opt.radiusDecrease;
end
% initialize
u      = u0;
lambda = opt.lambdaInit;
du = 0;
cap = @(x)max(0, min(1, x));
h   = initializeHistory(opt);
it   = 0;
if opt.plotEvolution
    fig = figure; 
end
accept = true;
radius = nan; 
while ~converged(it, h, opt)
    it   = it+1;
    uNew = cap(u+du);
    if opt.verbose
        [resNew, JNew] = f(uNew);
    else
        [~, resNew, JNew] = evalc('f(uNew)');
    end
    val = sum(resNew.^2);
    % update history for this step
    [h.val(it), h.lambda(it)] = deal(val, lambda);
    [h.u{it}, h.nIt(it)] = deal(uNew, h.nIt(it) + 1);
    if it > 1
        % actual to predicted improvement ratio 
        rho = -(val-h.val(it-1))/(dur'*(lambda*Dr*dur - gr));
        h.rho(it) = rho;
        accept = rho > 0;
        if rho < opt.ratioThresholds(1)
            lambda = lambda*badFac;
            radius = radius/badFac;
        elseif rho > opt.ratioThresholds(2)
            lambda = lambda/goodFac;
            radius = radius*goodFac;
        end
    end
    if it > 1 && ~accept
        % reset iteration, recompute update with increased lambda
        it = it-1;
    else
        % iteration accepted, compute new update
        [u, r, J] = deal(uNew, resNew, JNew);
        g = J'*r;                               % gradient
        % check norm of projected gradient
        pg = norm(u - cap(u-g));
        h.pg(it) = pg;
        if pg < opt.gradTol || it >= opt.maxIt +1
            continue;
        end
        isFree = ~(u==0 & g>0) & ~(u==1 & g<0); % free elems
        Jr = J(:, isFree);                      % reduced Jacobian
        gr = g(isFree);                         % reduced gradient
        
        JJr  = Jr'*Jr;
        if opt.scaledDamping 
            dr   = diag(JJr);
            mval = 1e-3*max(dr);
            dr(dr<mval) = mval;
            Dr   = diag(dr);                    % scaled diagonal
        else
            Dr = eye(nnz(isFree));              % identity
        end
    end
    [dur, lambda, radius] = computeUpdate(JJr, Dr, gr, lambda, radius, opt);
    du  = zeros(size(u));
    du(isFree) = dur;                                   % update
    if opt.plotEvolution
        fig = plotInfo(fig, it, h);
    end
    h.du(it) = norm(du);
    printInfo(h, it+~accept);
end 
[v, u, history] = handleOutput(it, h);
end

%--------------------------------------------------------------------------
function flag = converged(it, h, opt)
flag = false;
if it > 0
    dv = inf;
    if it > 1, dv = abs((h.val(it) - h.val(it-1))./h.val(it)); end
    flags = ...
        [it >= opt.maxIt+1, ...
        sum(h.nIt) >= opt.maxFunEvals, ...
        h.pg(it) < opt.gradTol, ...
        h.du(it) < opt.updateTol, ...
        h.val(it) < opt.resTolAbs, ...
        h.val(it)/h.val(1) < opt.resTolRel, ...
        dv < opt.resChangeTolRel];
        flag = any(flags);
end
if flag
    ll = 65;
    prnt = @(s)fprintf('| %s%s |\n', s, repmat(' ', [1, ll-length(s)-4]));
    fprintf('%s\n', repmat('-', [1, ll]));
    prnt('Optimization finnished:');
    switch find(flags, 1)
        case 1
            s = sprintf('Reached maximal number of iterations (%d)', it);
        case 2
            s = sprintf('Reached maximal number of function evaluations (%d)', sum(h.nIt(it)));
        case 3
            s = sprintf('Norm of projected gradient below tollerance (%7.2e < %7.2e)', h.pg(it), opt.gradTol);
        case 4 
            s = sprintf('Norm of update below tollerance (%7.2e < %7.2e)', h.du(it), opt.updateTol);
        case 5
            s = sprintf('Absolute mismatch below tollerance (%7.2e < %7.2e)', h.val(it), opt.resTolAbs);
        case 6
            s = sprintf('Relative mismatch below tollerance (%7.2e < %7.2e)', h.val(it)/h.val(1), opt.resTolRel);
        case 7
            s = sprintf('Relative mismatch change is below tollerance (%7.2e < %7.2e)', dv, opt.resChangeTolRel);
    end
    prnt(s)
    fprintf('%s\n', repmat('-', [1, ll]));
end
end

%--------------------------------------------------------------------------
function [du, lam, r] = computeUpdate(JJ, D, g, lam, r, opt)
if opt.lsqTol > eps
    lsqSolve = @(A,b)lsqminnorm(A, b, opt.lsqTol); 
else
    % backslash is probably OK and more efficient
    lsqSolve = @(A,b)mldivide(A,b);
end
if ~strcmp(opt.updateStrategy, 'TR')
    % simple standard approach
    lam = max(opt.lambdaMin, min(opt.lambdaMax, lam));
    du = -lsqSolve(JJ + lam*D, g);
else
    if ~isfinite(r)
        du  = -lsqSolve(JJ + lam*D, g);
        r   = norm(du);
    else
        % find lambda corresponding to r
        it  = 0;
        ndu = inf;
        lam0 = lam;
        while abs(ndu-r) > .1*r && it < 20
            % iteration taken from Fletcher's book (Practical Methods of
            % Optimization, Second Edition, p. 106)
            it = it +1;
            du  = -lsqSolve(JJ + lam*D, g);
            dut = -lsqSolve(JJ + lam*D, du);
            ndu  = norm(du);
            ndut = du'*dut/ndu;
            lam  = lam + (1 - ndu/r)*ndu/ndut;
            lam = max(opt.lambdaMin, min(opt.lambdaMax, lam));
        end
        if it == 20
            lam = lam0*opt.opt.lambdaIncrease;
            du  = -lsqSolve(JJ + lam*D, g);
            r   = norm(du);
        end
    end
end
end
    
%--------------------------------------------------------------------------
function [v, u, h] = handleOutput(it, h)
% check that last computed value is the best, otherwise remove last it
if it == 1
    warning('Optimization finnished after first function evaluation');
elseif h.val(it) > h.val(it-1)
    it = it-1;
end
v = h.val(it);
u = h.u{it};
fn = fieldnames(h);
for k = 1:numel(fn)
    h.(fn{k}) = h.(fn{k})(1:it);
end
end

%--------------------------------------------------------------------------
function h = initializeHistory(opt)
n = opt.maxIt+1;
h = struct('val', nan(1, n), 'u', {cell(1, n)}, 'lambda', nan(1,n), ...
           'rho', nan(1, n), 'nIt', zeros(1, n), 'pg',    nan(1,n), ...
           'du',  nan(1, n));
end
  
%--------------------------------------------------------------------------
function fig = plotInfo(fig, it, hst)
if ~ishandle(fig)
    fig = figure;
else
    % Avoid stealing focus if figure already exists
    set(0, 'CurrentFigure', fig);
end

xt = 0:(it-1);
xlim = [-.2, xt(end)+.5];
ch = abs(hst.val(2:end)-hst.val(1:end-1));
popt = {'o-', 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerFaceColor', [1 1 1]};
subplot(5,1,1)
semilogy(xt, abs(hst.val(1:it)), popt{:});
title('Objective');
set(gca, 'XLim', xlim)
subplot(5,1,2), semilogy(xt,hst.pg(1:it), popt{:}), title('Gradient norm');
set(gca, 'XLim', xlim)
subplot(5,1,3), semilogy(xt(2:end)-.5,ch(1:it-1), popt{:}), title('Objective change');
set(gca, 'XLim', xlim)
subplot(5,1,4), bar(xt,hst.nIt(1:it)), title('Local iterations');
set(gca, 'XLim', xlim)
subplot(5,1,5), col = [0.8500 0.3250 0.0980];
if numel(hst.u{it}) < 50
    bar(hst.u{it}, 'FaceColor', col);
else
    plot(hst.u{it}, '.-', 'Color', col);
end
title('Current scaled parameters');
set(gca, 'YLim', [0, 1])
drawnow
end    

%--------------------------------------------------------------------------
function printInfo(h, it)
fprintf('It: %2.1d | val: %4.3e | its: %3.1d | lambda: %4.3e | pgrad: %4.3e\n', ...
        it-1, h.val(it), h.nIt(it), h.lambda(it), h.pg(it));
end
