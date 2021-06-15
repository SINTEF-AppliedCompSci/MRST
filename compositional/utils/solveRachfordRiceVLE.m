function L = solveRachfordRiceVLE(L, K, z, varargin)
% Solve Rachford Rice equations to find liquid and vapor
% distribution.

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

    opt = struct('tolerance', 1e-12, ...
                 'min_z', 1e-10,...
                 'useAD', false, ...
                 'maxIterations', 1000);
    opt = merge_options(opt, varargin{:});
    if opt.useAD
        tmp1 = warning('query','MATLAB:nearlySingularMatrix');
        tmp2 = warning('query','MATLAB:singularMatrix');
        warning('off','MATLAB:nearlySingularMatrix')
        warning('off','MATLAB:singularMatrix')
    end
    % Ensure minimum value - equations degenerate when any component = 0.
    z = max(z, opt.min_z);
    nz = size(z, 1);
    % Allow special case - single set of K-values for all components
    if size(K, 1) < nz && size(K, 1) == 1
        K = repmat(K, size(z, 1), 1);
    end
    % Lowest bound for solution
    V_min = 1./(1 - max(K, [], 2));
    % Largest bound for solution
    V_max = 1./(1 - min(K, [], 2));
    if isempty(L)
        % Initial guess is the middle of the interval
        V = (V_min + V_max)/2;
    else
        V = 1 - L;
        bad_init = isnan(V);
        if any(bad_init)
            % NaN means that we do not know. Initialize to midpoint.
            V(bad_init) = (V_min(bad_init) + V_max(bad_init))/2;
        end
    end
    sz = [size(V, 1), size(K, 1), nz];
    assert(all(sz == max(sz)), ...
        'K, L (if not empty) and z inputs must all have equal number of rows.');
    n_V = numel(V);
    V_final = V;
    active = true(n_V, 1);
    for it = 1:opt.maxIterations
        % fprintf('Iteration %d: %d active\n', it, nnz(active));
        V0 = V;
        K_local = K(active, :);
        z_local = z(active, :);
        
        [dV, r] = getIncrement(V, K_local, z_local, opt);

        vNorm = abs(r);
        % Converged values do not need to be updated
        convResidual = vNorm < opt.tolerance;
        dV = ~convResidual.*dV;
        V = double(V) + dV;
        dLNorm = abs(V - V0);

        atBoundary = (dV > 0 & V0 == V_max) | (dV < 0 & V0 == V_min);
        convResidual(atBoundary) = true;
        conv = dLNorm < opt.tolerance | convResidual;

        bad = (V < V_min | V > V_max);
        if any(bad)
            % We are outside the bounds. Newton's method has failed. We
            % switch to a bisection, which is unconditionally convergent in
            % this internal (the objective function is monotone between the
            % singularities)
            V(bad) = bisection(V_max(bad), V_min(bad), K_local(bad, :),...
                               z_local(bad, :), opt);
            conv(bad) = true;
        end
        V = max(V, V_min);
        V = min(V, V_max);

        if it == opt.maxIterations
            warning('Reached max iterations in Rachford-Rice VLE calculation.')
            conv = conv | true;
        end
        update = false(n_V, 1);
        update(active) = conv;
        V_final(update) = V(conv);
        active(update) = false;
        V = V(~conv);
        V_max = V_max(~conv);
        V_min = V_min(~conv);

        if all(conv)
            break
        end
    end
    V = min(max(V_final, 0), 1);
    L = 1 - V;
    if opt.useAD
        warning(tmp1.state,'MATLAB:nearlySingularMatrix')
        warning(tmp2.state,'MATLAB:singularMatrix')
    end
end

function V = bisection(V_max0, V_min0, K, z, opt)
    fn = @(V) sum(((K - 1).*z)./(1 + bsxfun(@times, V, K - 1)), 2);
    left0 = fn(V_min0);
    right0 = fn(V_max0);
    [left, right, V_max, V_min] = deal(left0, right0, V_max0, V_min0);
    it = 1;
    tol = opt.tolerance;
    while it < 10000
        V = (V_max + V_min)/2;
        mid = fn(V);
        
        toLeft = mid.*left > 0;
        
        left(toLeft) = mid(toLeft);
        V_min(toLeft) = V(toLeft);
        right(~toLeft) = mid(~toLeft);
        V_max(~toLeft) = V(~toLeft);
        if norm(mid, inf) < tol || norm(V_min - V_max, inf) < tol
            break
        end
        it = it + 1;
    end
    % Treat case where endpoints are outside of the domain
    bad = sign(left) == sign(right);
    lrg = abs(left0(bad)) < abs(right0(bad));
    V(bad) = lrg.*V_min0(bad) + ~lrg.*V_max0(bad);
end


function [dV, r] = getIncrement(V, K, z, opt)
    if opt.useAD
        ncomp = size(K, 2);
        % Old AD version
        % V = initVariablesADI(V);
        V = initVariablesAD_diagonal(V);
        eq = 0;
        for i = 1:ncomp
            Ki = K(:, i);
            zi = z(:, i);
            eq = eq + ((Ki - 1).*zi)./(1 + V.*(Ki - 1));
        end                

        r = value(eq);
        if isa(eq, 'ADI')
            dV = -r./eq.jac{1}.diagonal;
            % J = -eq.jac{1};
            % dV = mldivide(J, r);
        else
            dV = 0*r;
        end
    else
        % Manual jacobian, faster
        if 1
            % Vectorized version
            mlt = @(x, y) bsxfun(@times, x, y);
            enum = sum((z.*(K - 1))./(1 + mlt(V, K-1)), 2);
            denom = sum((z.*(K - 1).^2)./(1 + mlt(V, K-1)).^2, 2);
        else
            ncomp = size(K, 2);
            % Partially vectorized version 
            enum = 0;
            denom = 0;
            for i = 1:ncomp
                Ki = K(:, i);
                zi = z(:, i);
                enum = enum + (zi.*(Ki - 1))./(1 + V.*(Ki-1));
                denom = denom + (zi.*(Ki - 1).^2)./(1 + V.*(Ki-1)).^2;
            end
        end
        dV = enum./denom;
        r = enum;
    end
end
