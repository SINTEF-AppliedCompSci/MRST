function c = findEnclosingCell(G, pt)
%Find cells with closest centroid (in Euclidian norm) in a 2D grid.
%
% SYNOPSIS:
%   c = findEnclosingCell(G, pt)
%
% PARAMETERS:
%   G  - Valid grid structure.  Must represent a two-dimensional geometry.
%
%   pt - Set of point coordinates, represented as an n-by-2 `double` array.
%
% RETURNS:
%   c  - Set of grid cells.  Specifically, c(i) for i=1:n, is the global
%        grid cell from 'G' for which ::
%
%           sum(bsxfun(@minus, G.cells.centroids, pt(i,:)) .^ 2, 2)
%
%        is minimised. If a point lies on the boundary between two cells,
%        the function returns the smallest index.
%
% SEE ALSO:
%   `pebi`.

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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


   assert (all(diff(G.faces.nodePos) == 2),                        ...
          ['Function ''findEnclosingCell'' is only supported in ', ...
           'two space dimensions.']);

   assert (size(pt, 2) == 2, ...
           'Nodal coordinates must be two-dimensinal.');

   cellno = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   sgn    = -1 + 2*(cellno==G.faces.neighbors(G.cells.faces(:,1), 1));
   edges = reshape(G.faces.nodes, 2, [])';

   % For each edge, check which side x is.
   x = G.nodes.coords(:,1);
   y = G.nodes.coords(:,2);

   a = [-diff(y(edges), 1, 2), diff(x(edges), 1, 2)];
   c = zeros([size(pt, 1), 1]);

   for k = 1 : size(pt, 1),
      b = bsxfun(@minus, pt(k,:), [x(edges(:,1)), y(edges(:,1))]);

      v = sum(a.*b, 2);
      V = ~ (sgn .* v(G.cells.faces(:,1)) < 0);

      i = accumarray(cellno, V, [G.cells.num, 1], @(x) all(x));
      ind = find(i);
      if isempty(ind)
         c(k) = 0;
      else
         c(k) = min(ind);
      end
   end
end
