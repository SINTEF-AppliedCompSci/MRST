function drawCellConnections3D(Gt, cells, varargin)
%Undocumented Utility Function

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

    x = Gt.cells.centroids(cells, 1);
    y = Gt.cells.centroids(cells, 2);
    z = Gt.cells.z(cells);
    smooth=@(x) x;
    x = smooth(x);
    y = smooth(y);
    z = smooth(z);
    plot3(x, y, z, varargin{:});
end
