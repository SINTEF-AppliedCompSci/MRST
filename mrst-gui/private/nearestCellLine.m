function [c f] = nearestCellLine(G, bf, pts)
% Find the cell nearest to a line defined by two points. Undocumented helper.

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

    N = numel(bf);
    x1 = repmat(pts(1,:), N, 1);
    x2 = repmat(pts(2,:), N, 1);
    x0 = G.faces.centroids(bf,:);
    % Find distance projected on line with a fudge factor for strange grids
    linedist = 10*sum(cross(x0-x1, x0-x2 ,2).^2, 2);
    ptdist = sum((x0-x1).^2, 2);
    dist = sqrt(ptdist - linedist);

    [val ind] = min(dist); %#ok backwards compatability
    f = bf(ind);
    c = G.faces.neighbors(f, :);
    c = c(c~=0);
end
