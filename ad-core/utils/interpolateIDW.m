function v = interpolateIDW(x, f, xq, order)
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

    % Number of sample points and dimension of space
    [ns, nx] = size(x);
    % Dimension of output space
    nf = size(f, 2);
    % Number of query points
    nq = size(xq, 1);
    
    assert(size(f, 1) == ns);
    assert(size(xq, 2) == nx);
    
    v = nan(nq, nf);
    for i = 1:nq
        d = bsxfun(@minus, x, xq(i, :));
        dist = max(sum(abs(d).^order, 2), 1e-20);
        wi = 1./dist;
        W = repmat(wi, 1, nf);
        for j = 1:nf
            disabled = isnan(f(:, j));
            W(disabled, j) = 0;
            f(disabled, j) = 0;
            s = sum(W(:, j));
            if s < 1e-18
                W(:, j) = repmat(1/ns, ns, 1);
            else
                W(:, j) = W(:, j)./s;
            end
        end
        v(i, :) = sum(W.*f, 1);
    end
end
