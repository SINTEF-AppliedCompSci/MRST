function [basis, X_bar] = createPod(values, frac_e, n_e)
%Undocumented Utility Function

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

    [X, X_bar] = snapshots(values);
    [V,S,W] = svd(X); % bruk econ?
    eig = diag(S).^2;
    tot_e = sum(eig);
    n = 0;
    energy = 0;
    while true
        n = n + 1;
        energy = energy + eig(n);
%         energy = eig(n);
        if energy/tot_e > frac_e || n == numel(eig) || n == n_e
            break
        end
    end

    fprintf('%d of %d eigenvalues included\n', n, numel(eig));
    basis = V(:, 1:n);
end
