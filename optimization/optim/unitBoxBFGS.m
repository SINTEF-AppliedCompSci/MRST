function [v, u, history] = unitBoxBFGS(u0, f, varargin)
% Iterative line search optimization using BFGS intended for scaled 
% problems, i.e., 0<=u<=1 and f~O(1) 
% 
% SYNOPSIS:
% [v, u, history] = unitBoxBFGS(u0, f, opt)

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
opt = struct(   'gradTol',             1e-3, ...
                'objChangeTol',        5e-4, ...
                'maxIt',               100,   ...
                'lineSearchMaxIt',     5,   ...
                'stepInit',            -1,   ...
                'wolfe1',              1e-3, ...
                'wolfe2',              0.9,  ...
                'safeguardFac',        1e-5, ...
                'stepIncreaseTol',     10,   ...
                'useBFGS',             true, ...
                'linEq',                 [], ...
                'linIneq',               [], ...
                'plotEvolution',       true);
opt  = merge_options(opt, varargin{:});
step = opt.stepInit;

if numel(u0) > 1e3
    warning('''unitBoxBFGS'' uses non-sparse representation of Hessian and constraints, and is no recommended for large (~>1000) number of controls')
end
% Perform initial evaluation of objective and gradient:
[v0, g0] = f(u0);
[v,u] = deal(v0,u0);
% If not provided, set initial step 
if step <= 0, step = 1; end
% Initial Hessian-approximation
Hi = -step*eye(numel(u0));
HiPrev = Hi;
% Setup constraint struct
c = getConstraints(u0, opt); 
% Setup struct for gathering optimization history
history = [];
% name|obj.val.|contr.|norm proj.grad.|ls-step|ls-its|ls-flag|hessian 
history = gatherInfo(history, v0, u0, nan, nan, nan, nan, Hi);
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
                r = 1/(du'*dg);
                V = eye(numel(u)) - r*dg*du';
                Hi = V'*Hi*V + r*(du*du');
            else
                fprintf('Hessian not updated during iteration %d.\n', it)
            end
        end
        % Update history
        history = gatherInfo(history, v, u, norm(pg,inf), lsinfo.step, ...
                             lsinfo.nits, lsinfo.flag, Hi);
    else
        history = gatherInfo(history, v, u, norm(pg,inf), 0, 1, 0, Hi);
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
% Equality constraints (always active)
if ~isempty(opt.linEq)
    sc = norm(opt.linEq.A);
    c.e.A = opt.linEq.A/sc;
    c.e.b = opt.linEq.b/sc;
else
    c.e.A = zeros(0, nu);
    c.e.b = [];
end
end
%--------------------------------------------------------------------------

function hst = gatherInfo(hst, val, u, pg, alpha, lsit, lsfl, hess)   
% obj.val | contr | norm proj.grad | ls-step | ls-its | ls-flag | hessian 
if isempty(hst)
    hst = struct('val', val, 'u', {{u}}, 'pg', pg, ...
                 'alpha', alpha, 'lsit', lsit, 'lsfl', lsfl, ...
                 'hess', {hess});
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
ch = [0, hst.val(2:end)-hst.val(1:end-1)];
subplot(5,1,1), plot(xt, hst.val, '.-','LineWidth', 2, 'MarkerSize', 20), title('Objective');
subplot(5,1,2), semilogy(xt,hst.pg,'.-','LineWidth', 2, 'MarkerSize', 20), title('Gradient norm');
subplot(5,1,3), semilogy(xt,ch,'.-','LineWidth', 2, 'MarkerSize', 20), title('Objective change');
subplot(5,1,4), bar(xt,hst.lsit), title('Line search iterations');
subplot(5,1,5), bar(hst.u{end}), title('Current scaled controls');
drawnow
end
%--------------------------------------------------------------------------

function [d, Hi, pg] = getSearchDirection(u0, g0, Hi, HiPrev, c)
% Find search-direaction which is the projection of Hi*g0 restricted to
% controls with non-active constraints. Check that direction is
% increasing, if not try HiPrev, if still not increasing, set Hi = I.

cnt = 1;
for k = 1:3
    if k==2
        Hi = HiPrev;
    elseif k==3
        Hi = -1;
    end
    % Check for active inequality constraints and project
    [na, na_prev] = deal(0, -inf);

    %[Q,~] = qr(c.e.A', 0);
    
    [Q, s] = svd(c.e.A', 0);
    s = diag(s);
    if ~isempty(s)
       Q = Q(:, s> sqrt(eps) * s(1));
    end

    d = - projQ(g0, Q, Hi);
    
    while na > na_prev
        ac = and(c.i.A*u0>=c.i.b-sqrt(eps), c.i.A*d >= -sqrt(eps));
    
        %[Q,~]  = qr([c.i.A(ac,:)', c.e.A'], 0);

        [Q, s] = svd([c.i.A(ac,:)', c.e.A'], 0);
        s = diag(s);
        if ~isempty(s)
           Q = Q(:, s > sqrt(eps)*s(1));
        end

        d  = - projQ(g0, Q, Hi);
        na_prev = na;
        na = nnz(ac);
    end
    pg = projQ(g0, Q);
    d = 0;
    done = false;
    gr = g0;
    while ~done
        dr      = - projQ(gr, Q, Hi);
        [ix, s] = findNextCons(c.i.A, c.i.b, u0+d, dr, ac);
        if ~isempty(ix);
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
    H = eye(numel(v));
end
    tmp = H*(v-Q*(Q'*v));
    w   = tmp - Q*(Q'*tmp);
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