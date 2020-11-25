function [G, cornercelltbl] = gridForPaperConvTest(Nx, gridType)
% We split the square (cube in 3d) in 4 (8 in 3d) parts. Indices of the cells of
% the (upper) north-west are returned.
    
% Planned Grid type:
% 1: Cartesian
% 2: Triangles by alternating bisection of triangles
% 3: Equilateral triangles
% 4: Triangles by uniform bisection (grid greated by Dolfin)

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}

    
    Nd = numel(Nx);
    
    switch gridType
      case 1
        G = computeGeometry(cartGrid(Nx, ones(1, Nd)));
        c = G.cells.centroids;
        cornercelltbl.cells = find(all(c > 0.5, 2));
        cornercelltbl = IndexArray(cornercelltbl);
      otherwise
        error('gridType not recognized');
    end
end
