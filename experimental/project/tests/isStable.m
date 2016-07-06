function stable = isStable(sol, g)

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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


tol = 1e-6;

stable = false;

magn = @(v)(sqrt(sum(v.^2,2)));

n = g.cells.normals(:,3)./magn(g.cells.normals);

h_z = n(:,1).*sol.h;

ix = (h_z < (g.cells.H-g.cells.H*tol)) & (h_z > g.cells.H*tol);

level = h_z + g.cells.z + abs(min(g.cells.z));

u = unique(level(ix));


if (max(u) < min(u) + max(u)*tol)
   stable = true;

end
