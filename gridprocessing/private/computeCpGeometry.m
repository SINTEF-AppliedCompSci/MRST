function G = computeCpGeometry(G)
%Compute Basic Geometric Properties Using Initial Corner-Point Coordinates
%
% SYNOPSIS:
%   G = computeCpGeometry(G)
%
% DESCRIPTION:
%   Computes cell and face centroids and cell extents using arithmetic
%   averages of initial corner points.  Sub-face geometry resulting from
%   faults is not considered here.
%
% PARAMETERS:
%   G - Grid structure.  Must provide corner vertices for all cells in an
%       array named `G.cells.cpnodes`.  Typically formed by construction
%       function `processGRDECL`.
%
% RETURNS:
%   G - Grid structure with added geometric information:
%         * cells.cpgeometry
%            - centroids: `G.cells.num` by `size(G.nodes.coords,2)` array
%              of approximate cell centre coordinates.
%
%            - extent: `G.cells.num` by `G.griddim` array of approximate
%              cell extents.  Derived from Euclidian distances between face
%              centroids on opposing sides.
%
%            - facecentroids: `size(G.cells.faces,1)` by
%              `size(G.nodes.coords,2)` array of approximate face centre
%              coordinates, relative to the cells.
%
% SEE ALSO:
%   `computeGeometry`, `processGRDECL`, `getFaceTransmissibility`.

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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

   assert (isfield(G.cells, 'cpnodes'), ...
          ['Corner-point geometry is not supported without explicit ', ...
           'mapping of a cell''s intial corner vertices']);

   [fc, cc] = derive_centroids(G);

   G.cells.cpgeometry = assign_cell_geometry(cellextent(G, fc), cc);
   G.cells.cpgeometry.facecentroids = assign_face_geometry(G, fc);
end

%--------------------------------------------------------------------------

function [fc, cc] = derive_centroids(G)
   % 3D array of size [ #cells, #coordinate directions, #cellnodes ]
   %
   % Cell index cycling faster than coordinate direction, coordinate
   % cycling faster than cell node ID.  Summing across dimension 3 thus
   % sums over the nodes in a node list.
   node_coords      = corner_points(G);
   make_face_centre = ...
      @(nodeID) sum(node_coords(:, :, nodeID), 3) ./ numel(nodeID);

   facenodes   = nominal_cell_faces(G);
   cf_per_cell = size(facenodes, 1);

   fc = zeros([G.cells.num * cf_per_cell, size(G.nodes.coords, 2)]);

   for k = 1 : cf_per_cell
      fc(k : cf_per_cell : end, :) = make_face_centre(facenodes(k, :));
   end

   % Put cell centre midway between the face centroids in the X direction.
   cc = (fc(1 : cf_per_cell : end, :) + fc(2 : cf_per_cell : end, :)) ./ 2;
end

%--------------------------------------------------------------------------

function cellfacecentroids = assign_face_geometry(G, fc)
   % Zero-based start index for cell's faces + facetag as index into fc.
   %
   % Explicitly relies on 'fc' being 2*d centroids-one for each of the 2*d
   % directional/nominal faces-per cell.

   cf_per_cell = cardinal_faces_per_cell(G);

   fIx = rldecode(0 : cf_per_cell : cf_per_cell * (G.cells.num - 1), ...
                  diff(G.cells.facePos), 2) .' ...
       + G.cells.faces(:, 2);

   % Broadcast corner-point face centroids to all pertinent faces.
   cellfacecentroids = fc(fIx, :);
end

%--------------------------------------------------------------------------

function geom = assign_cell_geometry(extent, cc)
   geom = struct('centroids', cc, 'extent', extent);
end

%--------------------------------------------------------------------------

function extent = cellextent(G, fc)
   distnorm    = create_distancevector_norm_function();
   cf_per_cell = cardinal_faces_per_cell(G);

   extent = zeros([G.cells.num, G.griddim]);

   for d = 1 : G.griddim
      % Extent in direction 'd' is Euclidian norm of distance between
      % centre points of opposing sides in Cardinal direction 'd'.
      extent(:, d) = distnorm(fc(2*d - 1 : cf_per_cell : end, :), ...
                              fc(2*d - 0 : cf_per_cell : end, :));
   end
end

%--------------------------------------------------------------------------

function node_coords = corner_points(G)
   % Size [#cells, #cellnodes, #coordinate directions]
   node_coords = reshape(G.nodes.coords(G.cells.cpnodes(:), :), ...
                         G.cells.num, size(G.cells.cpnodes, 2), []);

   % Swap dimensions 2 and 3 to create return value of size
   %
   %   [ #cells, #coordinate directions, #cellnodes ]
   %
   % Cell index cycling faster than coordinate direction, coordinate
   % cycling faster than cell node ID.
   node_coords = permute(node_coords, [1, 3, 2]);
end

%--------------------------------------------------------------------------

function facenodes = nominal_cell_faces(G)
   switch G.griddim
      case 1
         % Two faces, each constituted by a single vertex.
         facenodes = [ 1 ;
                       2 ];

      case 2
         % Four faces, each constituted by two vertices.
         facenodes = [ 1, 3 ;
                       2, 4 ;
                       1, 2 ;
                       3, 4 ];

      case 3
         % Six faces, each constituted by four vertices.
         facenodes = [ 1, 3, 5, 7 ;
                       2, 4, 6, 8 ;
                       1, 2, 5, 6 ;
                       3, 4, 7, 8 ;
                       1, 2, 3, 4 ;
                       5, 6, 7, 8 ];

      otherwise
         error('Griddim:Unsupported', ...
               'Unsupported grid dimension ''%d''', G.griddim);
   end
end

%--------------------------------------------------------------------------

function ncf = cardinal_faces_per_cell(G)
   ncf = 2 * G.griddim;
end

%--------------------------------------------------------------------------

function distvnorm = create_distancevector_norm_function()
   if exist('vecnorm', 'builtin')
      distvnorm = @(x1, x2) vecnorm(x2 - x1, 2, 2);
   else
      % Fall-back for 'vecnorm'; introduced in MATLAB 9.3.0 (R2017b)
      distvnorm = @(x1, x2) sqrt(sum((x2 - x1).^2, 2));
   end
end
