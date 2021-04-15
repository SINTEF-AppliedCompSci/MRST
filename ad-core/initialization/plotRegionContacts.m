function plotRegionContacts(G, region)
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

    coord = G.nodes.coords;
    xmin = min(coord(:, 1));
    xmax = max(coord(:, 1));
    ymin = min(coord(:, 2));
    ymax = max(coord(:, 2));
    X = [xmax, xmin, xmin, xmax];
    Y = [ymax, ymax, ymin, ymin];
    nc = numel(region.contacts);
    colors = jet(nc);


    for i = 1:nc
        c = region.contacts(i);
        Z = [c, c, c, c];
        patch(X, ...
              Y, ...
              Z, colors(i, :));
    end
end
