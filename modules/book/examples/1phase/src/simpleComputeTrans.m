function hT = simpleComputeTrans(G, rock)
%Compute transmissibilities.
%
% SYNOPSIS:
%   hT = simpleComputeTrans(G, rock)
%
% PARAMETERS:
%   G    - Grid structure as described by grid_structure.
%
%   rock - Rock data structure with valid field 'perm'.  The permeability
%          is assumed to be in measured in units of metres squared (m^2).
%          Use function 'darcy' to convert from darcies to m^2, e.g.,
%
%                 perm = convertFrom(perm, milli*darcy)
%
%          if the permeability is provided in units of millidarcies.
%
%          The field rock.perm may have ONE column for a scalar
%          permeability in each cell, TWO/THREE columns for a diagonal
%          permeability in each cell (in 2/3 D) and THREE/SIX columns for a
%          symmetric full tensor permeability.  In the latter case, each
%          cell gets the permeability tensor
%
%                 K_i = [ k1  k2 ]      in two space dimensions
%                       [ k2  k3 ]
%
%                 K_i = [ k1  k2  k3 ]  in three space dimensions
%                       [ k2  k4  k5 ]
%                       [ k3  k5  k6 ]
%
%
%
% RETURNS:
%   hT - half-transmissibilities for each local face of each grid cell in
%        the grid.  The number of half-transmissibilities equals the number
%        of rows in G.cells.faces.
%
% COMMENTS:
%   PLEASE NOTE: Face normals are assumed to have length equal to the
%   corresponding face areas.  This property is guaranteed by function
%   'computeGeometry'.
%
% SEE ALSO:
%   computeGeometry, simpleIncompTPFA, darcy, permTensor.

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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

   % Vectors from cell centroids to face centroids
   % To this end, we first need to determine the map from faces to cell
   % number so that the correct cell centroid is subtracted from each face
   % centroid.
   hf2cn = gridCellNo(G);
   C = G.faces.centroids(G.cells.faces(:,1), :) - G.cells.centroids(hf2cn,:);

   % Normal vectors
   % Face normals are assumed to have length equal to the corresponding
   % face areas. To get the correct sign, we look at the neighbouring
   % information that describes which cells share the face: if the current
   % cell number is in the first column, the face normal has a positive
   % sign. If not, it gets a negative sign
   sgn = 2*(hf2cn == G.faces.neighbors(G.cells.faces(:,1), 1)) - 1;
   N   = bsxfun(@times, sgn, G.faces.normals(G.cells.faces(:,1),:));
   clear sgn;
   
   % Extract permeability tensor.
   [K, i, j] = permTensor(rock, G.griddim);
   assert (size(K,1) == G.cells.num, ...
      ['Permeability must be defined in active cells only.\n', ...
      'Got %d tensors, expected %d (== number of cells).'],   ...
      size(K,1), G.cells.num);

   % Compute T = C'*K*N / C'*C. Loop-based to limit memory use.
   hT = zeros(size(hf2cn));
   for k=1:size(i,2),
      hT = hT + C(:,i(k)) .* K(hf2cn,k) .* N(:,j(k));
   end
   clear K i j hf2cn N;
   hT = hT./ sum(C.*C,2);
   clear C;

   is_neg = hT < 0;
   if any(is_neg),
      dispif(mrstVerbose, ...
            ['\nWarning:\n\t%d negative transmissibilities.\n\t', ...
             'Replaced by absolute values...\n'], sum(is_neg));

      hT(is_neg) = -hT(is_neg);
   end
   clear is_neg
end
