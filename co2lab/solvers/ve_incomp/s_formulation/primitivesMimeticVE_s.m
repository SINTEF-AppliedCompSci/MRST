function [areas, centVecs, outNormals] = primitivesMimeticVE_s(G, cf, cellno, sgn)
% Internal helper for topSurfaceGrid. Used to override mimetic primitives
% for computeMimeticIP to enable use of regular MRST solvers with the VE
% s-formulation. Intentionally left undocumented - see computeMimeticIP and
% computeMimeticIPVE instead.

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
nc = [G.nodes.coords G.nodes.z];
fc = [G.faces.centroids G.faces.z];
cc = [G.cells.centroids G.cells.z];
e_ij = nc(G.faces.nodes(1:2:end-1),:) - nc(G.faces.nodes(2:2:end),:);
areas        = sqrt(sum(e_ij.^2,2));
areas        = areas(cf);
centVecs_tmp = fc(cf,:) - cc(cellno,:);
surfBasis    = cross(rldecode(G.cells.normals, diff(G.cells.facePos)), ...
                     centVecs_tmp);
cfaceNormals = zeros(size(surfBasis,1),2);
centVecs     = zeros(size(surfBasis,1),2);

for i=1:G.cells.num
   cfind = G.cells.facePos(i):G.cells.facePos(i+1)-1;
   C3D   = centVecs_tmp(cfind,:);
   ee3D  = e_ij(cf(cfind),:);
   U     = orth(surfBasis(cfind,:)')';
   if (dot(cross(U(1,:),U(2,:)), G.cells.normals(i,:)) < 0)
      U = U([2 1],:);
      assert(dot(cross(U(1,:),U(2,:)),G.cells.normals(i,:))>0);
   end
   C2D  = U*C3D';
   ee2D = (U*ee3D')';
   N2D  = [-ee2D(:,2),ee2D(:,1)];
   centVecs(cfind,:) = C2D';
   cfaceNormals(cfind,:) = N2D;
end

outNormals = bsxfun(@times, cfaceNormals, sgn);
end
