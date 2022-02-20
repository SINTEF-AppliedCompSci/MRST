function N = findStreamlineNeighborshipCart(G)
% Build (n x 2*d) -array of neighbors for each cell in (Cartesian) grid G.

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

   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   col    = 1 +   (cellNo == G.faces.neighbors(G.cells.faces(:,1), 1));
   c      = G.faces.neighbors(double(G.cells.faces(:,1)) + G.faces.num* (col-1));
   N      = accumarray([cellNo, G.cells.faces(:,2)], c);
end
