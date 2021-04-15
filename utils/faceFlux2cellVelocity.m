function v = faceFlux2cellVelocity(G, faceFlux)
%Transform face-based flux field to one constant velocity per cell.
%
% SYNOPSIS:
%   veclocity = faceFlux2cellVelocity(G, faceFlux)
%
% PARAMETERS:
%   G        - Grid structure.
%
%   faceFlux - Vector of fluxes corresponding to face ordering.
%
% RETURNS:
%   velocity - G.cells.num-by-d matrix of cell velocities.
%
% SEE ALSO:
%   `cellFlux2faceFlux`, `faceFlux2cellFlux`.

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


   % Compute constant part of velocity field in each cell

   if ~isfield(G.cells, 'centroids'),
      G = computeGeometry(G);
   end

   [cellNo, cellFaces] = getCellNoFaces(G);

   C  = G.faces.centroids(cellFaces, :) - ...
        G.cells.centroids(cellNo   , :);

   cf = faceFlux2cellFlux(G, faceFlux);

   v = sparse(cellNo, 1 : numel(cellNo), cf) * C;
   v = bsxfun(@rdivide, v, G.cells.volumes);
end
