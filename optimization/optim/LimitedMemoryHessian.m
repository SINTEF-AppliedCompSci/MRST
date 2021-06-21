classdef LimitedMemoryHessian
% LimitedMemoryHessian: simple implementation of limited memory
% approxiamtion of invere Hessian (for use in L-BFGS)
% see: unitBoxBFGS (option limitedMemory)
properties
    initScale       % scaling for H0 (defualt 1)
    initStrategy    % 'static' : H0k = H0 (default)
                    % 'dynamic': H0k = gamma_k*I  
    nullspace = []; % nullspace of active subspace 
end
properties (SetAccess = private)
    m               % number of vectors stored (default 5)
    S        = [];  % m recent control/param diffs
    Y        = [];  % m recent gradient diffs 
    itCount  = 0;   % update iteration count
    sign            % sign of Hessian(default 1)
end

methods 
    function H = LimitedMemoryHessian(varargin)
        opt = struct('initScale',           1, ...
                     'initStrategy', 'static', ...
                     'm',                   5, ...
                     'sign',                1);
        opt = merge_options(opt, varargin{:});
        fn = fieldnames(opt);
        for k = 1:numel(fn)
            H.(fn{k}) = opt.(fn{k});
        end
    end
    % ---------------------------------------------------------------------
    
    function H = update(H, s ,y)
        % store new vectors
        H.itCount = H.itCount +1;
        pix = ':';
        if H.itCount > H.m
            pix = 2:H.m;
        end
        H.S = [H.S(:, pix), s];
        H.Y = [H.Y(:, pix), y];
        if H.itCount == 1 % limit m to number of vars
            if H.m > numel(s)
                H.m = numel(s);
                dispif(mrstVerbose, 'Resetting ''m'' to number of parameters (%d)\n', H.m);
            end
        end 
    end
    
    % ---------------------------------------------------------------------
    function H = reset(H)
        [H.S, H.Y] = deal([]);
        H.itCount    = 0;
        H.nullspace   = [];
    end
    
    % ---------------------------------------------------------------------
    function r = mtimes(H, v)
        % apply mulitplication H*v 
        if isa(v, 'LimitedMemoryHessian')
            error('Only right-multiplciation is supported');
        else
            if H.itCount == 0
                r = (H.sign*H.initScale)*v;
                if ~isempty(H.nullspace)
                    r = r - H.nullspace*(H.nullspace'*r);
                end
            else
                assert(all(size(v)==[size(H.S, 1), 1]), ...
                    'Dimension mismatch');
                if isempty(H.nullspace) % do standard L-BFGS
                    nVec  = size(H.S, 2);
                    [rho, alpha] = deal(nan(1, nVec));
                    for k = nVec:-1:1
                        rho(k)   = 1/(H.S(:,k)'*H.Y(:,k));
                        alpha(k) = rho(k)*(H.S(:,k)'*v);
                        v        = v - alpha(k)*H.Y(:,k);
                    end
                    r = H.applyInitial(v);
                    for k = 1:nVec
                        beta = rho(k)*(H.Y(:,k)'*r);
                        r    = r + (alpha(k)-beta)*H.S(:,k);
                    end
                else % do subspace version
                    th = H.applyInitial(1);
                    r  = subspaceProd(H.S, H.Y, th, H.nullspace, v);
                end
            end
        end
    end
    
    % ---------------------------------------------------------------------
    function H = uminus(H)
        H.sign = -H.sign;
    end 
    
    % ---------------------------------------------------------------------
    function r = applyInitial(H, v)
        if strcmp(H.initStrategy, 'static')
            r = (H.sign*H.initScale)*v;
        elseif strcmp(H.initStrategy, 'dynamic')
            [s,y] = deal(H.S(:,end), H.Y(:,end));
            r = ((s'*y)/(y'*y))*v;
        else
            error('Unknown strategy: %s', H.initStrategy);
        end
    end
    
    % ---------------------------------------------------------------------
    function H = setNullspace(H, Q)
        if nargin < 2
            H.nullspace = [];
        else
            if H.itCount > 0
                assert(size(Q,1)==size(H.S, 1), 'Dimension mismatch');
                assert(size(Q,1)>=size(Q,2), ...
                    'Number of columns in nullspace matrix exceeds number of rows');
            end
            H.nullspace = Q;
        end
    end
    
    % ---------------------------------------------------------------------
    function M = full(H)
        % Create full matrix, only intended for debugging purposes
        if H.itCount == 0
            M = H*1;
        else
            n = size(H.S, 1);
            M = zeros(n);
            for k = 1:n
                ek = zeros(n,1);
                ek(k)  = 1;
                M(:,k) = H*ek;
            end
        end
    end
end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function r = subspaceProd(S, Y, th, Q, v)
% Perform product restricted to nullspace of Q-colunms, i.e.,
% r = Z * (Z'*Hessian*Z)^-1 * Z'*v
% with Z = null(Q')
% Formula from: 
% Byrd, R. H., Lu, P., Nocedal, J., & Zhu, C. (1995). 
% A limited memory algorithm for bound constrained optimization. 
% SIAM Journal on scientific computing, 16(5), 1190-1208.
%
% S  : control diffs
% Y  : gradient diffs
% th : inverse Hessian scaling
% Q  : nullspace of active subspace
% v  : vector
if size(Q,1) <= size(Q,2)
    % null(Q')=[], return zero;
    r = 0*v;
else
    n   = size(S,2);
    W   = [Y, S/th];
    tmp = S'*Y;
    L   = tril(tmp, -1);
    D   = diag(diag(tmp));
    M   = [-D, L'; L, (S'*S)/th];
    % projection onto active subspace (u -> Z*Z'*u)
    projSub = @(u)u-Q*(Q'*u);
    r   = M\(W'*projSub(v));
    r   = (eye(2*n) - th*(M\(W'*projSub(W))))\r;
    r   = projSub(th*v + th^2*(W*r));
end
end

%{
% For reference: 
% should be equal to subspaceProd, but at the additional cost of doing null(Q') 
Z = null(Q');
n   = size(S,2);
W   = [Y, S/th];
tmp = S'*Y;
L   = tril(tmp, -1);
D   = diag(diag(tmp));
M   = [-D, L'; L, (S'*S)/th];
WtZ = W'*Z;
Ztv  = Z'*v; 
r   = M\(WtZ*Ztv);
r   = (eye(2*n)-th*(M\(WtZ*WtZ')))\r;
r   = th*Ztv + (th^2)*WtZ'*r;
r   = Z*r;
%}

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
