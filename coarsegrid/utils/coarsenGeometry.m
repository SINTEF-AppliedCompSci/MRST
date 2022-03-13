function cg = coarsenGeometry(cg)
% Add geometry (centroids, face normals, areas, ...) to a coarse grid
%
% SYNOPSIS:
%   cg = coarsenGeometry(cg)
%
% PARAMETERS:
%   cg    - Coarse grid including parent (fine) grid. (New definition).
%
% RETURNS:
%   cg    - Coarse grid with accumulated areas or volumes, area- or
%           volume-weighted centroids and accumulated area-weighted
%           normals accounting for sign convention in coarse and fine grid.
%
% NOTE:
%   The coarse geometry is derived from the geometry on the underlying fine
%   grid. For this reason, the underlying fine grid should be the output
%   from a call to `computeGeometry(G)`.
%
% SEE ALSO:
%   `generateCoarseGrid`, `computeGeometry`

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


   assert (isfield(cg, 'parent'), ...
          ['Huh!? Field ''parent''missing in coarse grid.  ', ...
           'Did you really supply a coarse grid?']);

   assert (isfield(cg.parent.cells, 'volumes'), ...
          ['Geometric fields missing in fine grid.  ', ...
           'Did you forget to call computeGeometry on fine grid?']);

   % Cells
   nc = cg.parent.cells.num;
   c  = sparse(cg.partition, 1 : nc, cg.parent.cells.volumes) ...
        *                                                     ...
        [cg.parent.cells.centroids, ones([nc, 1])];

   cg.cells.volumes   = c(:, end);
   cg.cells.centroids = bsxfun(@rdivide, c(:, 1:end-1), cg.cells.volumes);

   clear c

   % Faces
   faceno = rldecode(1:cg.faces.num, diff(cg.faces.connPos), 2)';
   sgn    = fineToCoarseSign(cg);
   for i = 1:size(cg.parent.faces.centroids, 2),
      cg.faces.centroids(:,i) = accumarray(faceno, ...
         cg.parent.faces.centroids(cg.faces.fconn, i).*...
         cg.parent.faces.areas(cg.faces.fconn));

      cg.faces.normals(:,i) = accumarray(faceno, ...
         cg.parent.faces.normals(cg.faces.fconn, i).*sgn);
   end
   totalareas         = accumarray(faceno, cg.parent.faces.areas(cg.faces.fconn));
   cg.faces.centroids = bsxfun(@rdivide, cg.faces.centroids, totalareas);
   cg.faces.areas     = sqrt(sum(cg.faces.normals.^2,2));
   
end
