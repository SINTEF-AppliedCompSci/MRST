function L = solveRachfordRiceVLE(L, K, z, varargin)
% Solve Rachford Rice equations to find liquid and vapor
% distribution.
    opt = struct('tolerance', 1e-12, ...
                 'min_z', 1e-10,...
                 'maxIterations', 100);
    opt = merge_options(opt, varargin{:});
    tmp1 = warning('query','MATLAB:nearlySingularMatrix');
    tmp2 = warning('query','MATLAB:singularMatrix');
    warning('off','MATLAB:nearlySingularMatrix')
    warning('off','MATLAB:singularMatrix')

    maxit = 100;
    K(~isfinite(K)) = 1e6;
    
    singularities = 1./(1 - K);
    singularities = min(max(singularities, 0), 1);
    
    L_min = min(singularities, [], 2);
    L_max = max(singularities, [], 2);
    if isempty(L)
        L = (L_min + L_max)/2;
    end
    n_L = numel(L);
    L_final = L;
    active = true(n_L, 1);
    for it = 1:opt.maxIterations
        % fprintf('Iteration %d: %d active\n', it, nnz(active));
        L0 = L;
        K_local = K(active, :);
        z_local = z(active, :);
        
        [dL, r] = getIncrement(L, K_local, z_local, opt);

        vNorm = abs(r);
        % Converged values do not need to be updated
        convResidual = vNorm < opt.tolerance;
        dL = ~convResidual.*dL;
        % dL = sign(dL).*min(abs(dL), 0.2);
        L = double(L) + dL;
        dLNorm = abs(L - L0);

        atBoundary = (dL > 0 & L0 == L_max) | (dL < 0 & L0 == L_min);
        convResidual(atBoundary) = true;
        conv = dLNorm < opt.tolerance | convResidual;

        bad = (L < L_min | L > L_max);
        if any(bad)
            % Perform bisection
            L(bad) = bisection(L_max(bad), L_min(bad), K_local(bad, :), z_local(bad, :), opt);
            conv(bad) = true;
        end
        L = max(L, L_min);
        L = min(L, L_max);

        if it == maxit
            disp('Reached max iterations in Rachford rice')
            conv = conv | true;
        end
        update = false(n_L, 1);
        update(active) = conv;


        L_final(update) = L(conv);
        active(update) = false;
        L = L(~conv);
        L_max = L_max(~conv);
        L_min = L_min(~conv);

        if all(conv)
            break
        end
    end
    L = L_final;
    warning(tmp1.state,'MATLAB:nearlySingularMatrix')
    warning(tmp2.state,'MATLAB:singularMatrix')
end

function L = bisection(L_max0, L_min0, K, z, opt)
    fn = @(L) sum(((K - 1).*z)./(1 + (1-L).*(K - 1)).*(z>opt.min_z), 2);
    left0 = fn(L_min0);
    right0 = fn(L_max0);
    [left, right, L_max, L_min] = deal(left0, right0, L_max0, L_min0);
    it = 1;
    tol = opt.tolerance;
    while it < 10000
        L = (L_max + L_min)/2;
        mid = fn(L);
        
        toLeft = mid.*left > 0;
        
        left(toLeft) = mid(toLeft);
        L_min(toLeft) = L(toLeft);
        right(~toLeft) = mid(~toLeft);
        L_max(~toLeft) = L(~toLeft);
        if norm(mid, inf) < tol || norm(L_min - L_max, inf) < tol
            break
        end
        it = it + 1;
    end
    % Treat case where endpoints are outside of the domain
    bad = sign(left) == sign(right);
    lrg = abs(left0(bad)) < abs(right0(bad));
    L(bad) = lrg.*L_min0(bad) + ~lrg.*L_max0(bad);
end


function [dL, r] = getIncrement(L, K, z, opt)
    ncomp = size(K, 2);
    if 0
        % Old AD version
        L = initVariablesADI(L);
        eq = 0;
        for i = 1:ncomp
            Ki = K(:, i);
            zi = z(:, i);
            v = ((Ki - 1).*zi)./(1 + (1-L).*(Ki - 1));
            present = zi > opt.min_z;
            eq = eq + v.*present;
        end                

        r = double(eq);
        if isa(eq, 'ADI')
            J = -eq.jac{1};
            dL = mldivide(J, r);
        else
            dL = 0*r;
        end
    else
        % Manual jacobian, faster
        enum = 0;
        denom = 0;
        V = 1 - L;
        for i = 1:ncomp
            Ki = K(:, i);
            zi = z(:, i);
            present = zi > opt.min_z;

            enum = enum + present.*(zi.*(Ki - 1))./(1 + V.*(Ki-1));
            denom = denom + present.*(zi.*(Ki - 1).^2)./(1 + V.*(Ki-1)).^2;
        end
        dL = - enum./denom;
        r = enum;
    end
end