function p_f = interpolatedFacePressure(G, p, useDist)
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

    if nargin == 2
        useDist = true;
    end
    p_set = [0; p];
    n = G.faces.neighbors + 1;
    if useDist
        c = [nan(1, G.griddim); G.cells.centroids];
        
        d_f = zeros(G.faces.num, 1);
        for i = 1:2
            d_f(:, i) = 1./sqrt(sum((c(n(:, i), :) - G.faces.centroids).^2, 2));
        end
        d_f(isnan(d_f)) = 0;
        p_f = bsxfun(@rdivide, sum(p_set(n).*d_f, 2), sum(d_f, 2));
    else
        p_f = sum(p_set(n), 2)./sum(G.faces.neighbors > 0, 2);
    end
end
