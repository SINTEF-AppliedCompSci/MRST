classdef LimitedMemoryHessian
    % LimitedMemoryHessian: simple implementation of limited memory
    % approxiamtion of invere Hessian (for use in L-BFGS)
    % see: unitBoxBFGS (option limitedMemory)
properties
    initScale     % scaling for H0 (defualt 1)
    initStrategy  % 'static' : H0k = H0 (default)
                  % 'dynamic': H0k = gamma_k*I  
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
    end
    % ---------------------------------------------------------------------
    function H = reset(H)
        [H.S, H.Y] = deal([]);
        H.itCount    = 0;
    end
    % ---------------------------------------------------------------------
    function r = mtimes(H, v)
        % apply mulitplication H*v 
        if isa(v, 'LimitedMemoryHessian')
            error('Only right-multiplciation is supported');
        else
            if H.itCount == 0
                r = (H.sign*H.initScale)*v;
            else
                assert(all(size(v)==[size(H.S, 1), 1]), ...
                    'Dimension mismatch');
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
end
end
