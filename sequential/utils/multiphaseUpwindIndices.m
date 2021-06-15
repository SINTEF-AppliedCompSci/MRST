function [up, theta, r] = multiphaseUpwindIndices(G, vT, T, K, upstr)
    % Implementation of algorithm from "UPSTREAM DIFFERENCING FOR
    % MULTIPHASE FLOW IN RESERVOIR SIMULATIONS" by Brenier & Jaffre
    nPh = numel(G);
    nF  = numel(T);
    r = zeros(nF, 1);
    theta = repmat(vT, 1, nPh);
    if nF == 0
        up = true(nF, nPh);
        return
    end
    G = value(G);
    K = value(K);
    % Sort phases for each interface by potential
    [G, sortInd] = sort(G, 2);
    % Get the linear index for this sorting
    subs = sub2ind([nF, nPh], repmat((1:nF)', 1, nPh), sortInd);
    flag = true(nF, 1);
    left_upwind = upstr(flag, K);
    right_upwind = upstr(~flag, K);
    left_upwind = left_upwind(subs);
    right_upwind = right_upwind(subs);
    for l = 1:nPh
        for j = 1:nPh
            if j == l
                % Potential difference between a phase and itself is zero
                continue
            end
            if j > l
                % Phase has larger potential, take upwind
                kj = left_upwind(:, j);
            else
                % Phase has lower potential, take downwind
                kj = right_upwind(:, j);
            end
            theta(:, l) = theta(:, l) + T.*(G(:, l) - G(:, j)).*kj;
        end
    end
    % r is zero if all thetas are positive, otherwise it is the last phase
    % with negative value
    for l = 1:nPh
        r(theta(:, l) <= 0) = l;
    end
    up = false(nF, nPh);
    for l = 1:nPh
        up(:, l) = l > r;
    end
    % Permute back to original ordering
    up(subs) = up;
end

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
