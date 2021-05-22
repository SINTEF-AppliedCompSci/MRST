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
%       'maximize' : Boolean option (false will minimize objective)
%                    Default: true
%       'stepInit' : Initial step gradient scaling. If not provided, the 
%                    following scaling will be used:
%                    if 0.1 < |initial objective| < 10
%                        stepInit = 1
%                    otherwise
%                        stepInit = 0.1/max(|initial gradient|)
%   Stopping criteria options
%       'gradTol'       : Absolute tollerance of inf-norm of projected gradient. 
%                         Default: 1e-3
%       'objChangeTol'  : Absolute objective update tollerance. Default: 5e-4
%       'maxIt'         : Maximal number of iterations. Default: 25
%   Line-search options
%       'lineSearchMaxIt' : Maximal number of line-search iterations. Default: 5
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
%       'linEq'    : Linear equality constraints given as structure with
%                    fields 'A' and 'b' to enforce A*u=b
%       'linIneq'  : Linear inequality constraints *in addition to* the default 
%                    box constraints (0 <= u <= 1). Given as structure with
%                    fields 'A' and 'b' to enforce A*u<=b
%   Plotting and output options
%       'plotEvolution' : Plot progess of optimization in figure. 
%                         Default: true
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
                'gradTol',             1e-3, ...
                'objChangeTol',        5e-4, ...
                'maxIt',               100,   ...
                'lineSearchMaxIt',     5,   ...       
                'wolfe1',              1e-3, ...
                'wolfe2',              0.9,  ...
                'safeguardFac',        1e-5, ...
                'stepIncreaseTol',     10,   ...
                'useBFGS',             true, ...
                'limitedMemory',       false, ...
                'lbfgsNum',            5, ...
                'lbfgsStrategy',       'static', ...
                'linEq',               [], ...
                'linIneq',             [], ...
                'plotEvolution',       true, ...
                'outputHessian',       false);
opt  = merge_options(opt, varargin{:});
objSign = 1;
if ~opt.maximize
    f = @(u)fNegative(u, f);
    objSign = -1;
end

if numel(u0) > 5e2
    warning('Switching to L-BFGS with default options due to large number of controls')
    opt.limitedMemory = true;
end
% Perform initial evaluation of objective and gradient:
[v0,g0] = f(u0);
[v ,u ] = deal(v0,u0);
% If not provided, set initial step 
step = opt.stepInit;
if isnan(step)
    if v0>.1 && v0 < 10
        step = 1;
    else
        step = 0.1/max(abs(g0));
    end
end
% Initial Hessian-approximation
if ~opt.limitedMemory
    Hi = step*eye(numel(u0));
else
    Hi = LimitedMemoryHessian('initScale', step, ...
                              'm', opt.lbfgsNum, ...
                              'initStrategy', opt.lbfgsStrategy, ...
                              'sign', 1);
end
HiPrev = Hi;
% Setup constraint struct
c = getConstraints(u0, opt); 
% Setup struct for gathering optimization history
history = [];
% name|obj.val.|contr.|norm proj.grad.|ls-step|ls-its|ls-flag|hessian 
history = gatherInfo(history, objSign*v0, u0, norm(g0), nan, nan, nan, Hi);
if opt.plotEvolution
    plotInfo(10, history)
end
    
it = 0;
success = false;
g = g0;
while ~success
    it = it+1;
    [d, Hi, pg] = getSearchDirection(u0, g0, Hi, HiPrev, c);
    if isempty(d), break; end
    if ~(norm(pg,inf) < opt.gradTol)
        % Perform line-search
        [u, v, g, lsinfo] = lineSearch(u0, v0, g0, d, f, c, opt);
        
        % Update Hessian approximation
        if opt.useBFGS && lsinfo.flag == 1
            [du, dg] = deal(u-u0, g-g0);
            if abs(du'*dg) > sqrt(eps)*norm(du)*norm(dg)
                HiPrev = Hi;
                if isa(Hi, 'LimitedMemoryHessian')
                    Hi = Hi.update(du, dg);
                else
                    r = 1/(du'*dg);
                    V = eye(numel(u)) - r*dg*du';
                    Hi = V'*Hi*V + r*(du*du')
                end
            else
                fprintf('Hessian not updated during iteration %d.\n', it)
            end
        end
        % Update history
        history = gatherInfo(history, objSign*v, u, norm(pg,inf), lsinfo.step, ...
                             lsinfo.nits, lsinfo.flag, Hi);
    else
        history = gatherInfo(history, objSign*v, u, norm(pg,inf), 0, 1, 0, Hi);
    end
    
    %Check stopping criteria
    success = (it >= opt.maxIt) || (norm(pg,inf) < opt.gradTol) || ...
              (abs(v-v0) < opt.objChangeTol);

    [u0, v0, g0] = deal(u, v, g);
    if opt.plotEvolution
        plotInfo(10, history)
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
[Q, s] = svd(c.i.A', 0);
s      = diag(s);
c.i.Q  = Q(:, s> sqrt(eps) * s(1));
c.i.isActive = false(size(c.i.b));
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

function hst = gatherInfo(hst, val, u, pg, alpha, lsit, lsfl, hess)   
% obj.val | contr | norm proj.grad | ls-step | ls-its | ls-flag | hessian 
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
  
function [] = plotInfo(fig, hst)
if ~ishandle(fig)
    figure(fig)
else
    % Avoid stealing focus if figure already exists
    set(0, 'CurrentFigure', fig);
end
xt = 0:(numel(hst.val)-1);
xlim = [-.2, xt(end)+.5];
ch = [0, abs(hst.val(2:end)-hst.val(1:end-1))];
subplot(5,1,1), plot(xt, hst.val, '.-','LineWidth', 2, 'MarkerSize', 20), title('Objective');
set(gca, 'XLim', xlim)
subplot(5,1,2), semilogy(xt,hst.pg,'.-','LineWidth', 2, 'MarkerSize', 20), title('Gradient norm');
set(gca, 'XLim', xlim)
subplot(5,1,3), semilogy(xt,ch,'.-','LineWidth', 2, 'MarkerSize', 20), title('Objective change');
set(gca, 'XLim', xlim)
subplot(5,1,4), bar(xt,hst.lsit), title('Line search iterations');
set(gca, 'XLim', xlim)
subplot(5,1,5), bar(hst.u{end}, 'FaceColor', 'g'), title('Current scaled controls');
set(gca, 'YLim', [0, 1])
drawnow
end
%--------------------------------------------------------------------------

function [d, Hi, pg] = getSearchDirection(u0, g0, Hi, HiPrev, c)
% Find search-direaction which is the projection of Hi*g0 restricted to
% controls with non-active constraints. Check that direction is
% increasing, if not try HiPrev, if still not increasing, set Hi = I.
sgn = 1;
cnt = 1;
for k = 1:3
    if k==2
        Hi = HiPrev;
    elseif k==3
        if isa(Hi, 'LimitedMemoryHessian')
            Hi = Hi.reset();
        else
            Hi = -1;
        end
    end
    % Check for active inequality constraints and project
    [na, na_prev] = deal(0, -inf);

    [Q, s] = svd(c.e.A', 0);
    s = diag(s);
    if ~isempty(s)
       Q = Q(:, s> sqrt(eps) * s(1));
    end

    d = sgn*projQ(g0, Q, Hi)
    
    while na > na_prev
        ac = and(c.i.A*u0>=c.i.b-sqrt(eps), c.i.A*d >= -sqrt(eps));
    
        %[Q,~]  = qr([c.i.A(ac,:)', c.e.A'], 0);

        [Q, s] = svd([c.i.A(ac,:)', c.e.A'], 0);
        s = diag(s);
        if ~isempty(s)
           Q = Q(:, s > sqrt(eps)*s(1));
        end

        d  = sgn*projQ(g0, Q, Hi);
        na_prev = na;
        na = nnz(ac);
    end
    pg = sgn*projQ(g0, Q);
    d = 0;
    done = false;
    gr = g0;
    while ~done
        dr      = sgn*projQ(gr, Q, Hi);
        [ix, s] = findNextCons(c.i.A, c.i.b, u0+d, dr, ac);
        if ~isempty(ix)
            ac(ix) = true;
                d  = d + s*dr;
                gr = (1-s)*gr;
                Q  = expandQ(Q, c.i.A(ix,:)');
        else
            d = d+dr;
            done = true;
        end
    end
    if d'*g0 >= 0 % increasing search direction found, exit
        break;
    else
        cnt = cnt+1;
        if cnt == 2
            fprintf('Non-increasing search direction, using previous Hessian approximation.\n')
        else
            fprintf('Non-increasing search direction, resetting Hessian to identity.\n')
        end
    end
end
if d'*g0 < 0
    fprintf('All controls on constraints or function decreasing along gradient\n')
    d = [];
end

if max(abs(d)) < sqrt(eps)
   fprintf('Gradient too small to allow computation ot steplength\n');
   d = [];
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

function [ix, s] = findNextCons(A, b, u, d, ac)
s = (b-A*u)./(A*d);
s(ac)  = inf;
s(s<eps) = inf;
[s, ix] = min(s);
if s >= 1
    ix = [];
end
end
%--------------------------------------------------------------------------

function Q = expandQ(Q, v)
v = v - Q*(Q'*v);
if norm(v) > 100*eps
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

