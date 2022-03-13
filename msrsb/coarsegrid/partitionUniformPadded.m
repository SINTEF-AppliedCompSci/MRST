function p = partitionUniformPadded(G, dims)
% Equivialent to partitionUI from the coarsegrid module, with a small
% half-block at each end of the domain.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
    ijk = gridLogicalIndices(G);
    ijk = [ ijk{:} ];
    M = max(ijk) - min(ijk) + 1;
    coarseDim = dims - 1;
    spacings = cell(G.griddim, 1);
    for d = 1:G.griddim
       v = ijk(:, d);
       B = coarseDim(d);
       if B == 0
           spacings{d} = M(d);
       else
           cts = lbLinDist((min(v) - 1):max(v) - 1, M(d), B)';
           [~, s] = rlencode(cts);

           divEl = s(1);
           e1 = floor(divEl/2);
           e2 = divEl - e1;
           spacings{d} = [e1; s(2:end); e2];
       end
    end
    p = partitionTensor(G, spacings{:});
end

function f = lbLinDist(f, M, B)
    % Stolen from partitionUI
    L = floor(M ./ B);  % Tentative number of cells per coarse block.
    R = mod(M, B);      % Additional cells not previously accounted for.
    f = max(floor(f ./ (L + 1)), floor((f - R) ./ L));
end