function [H, gcells, gfaces, gnodes] = extractSubgrid(G, c)
%Construct valid grid definition from subset of existing grid cells.
%
% SYNOPSIS:
%    subG              = extractSubgrid(G, cells)
%   [subG, gc, gf, gn] = extractSubgrid(G, cells)
%
% PARAMETERS:
%   G     - Valid grid definition.
%
%   cells - Subset of existing grid cells for which the subgrid data
%           structure should be constructed.  Must be an array of numeric
%           cell indices or a logical mask into the range 1:G.cells.num.
%           Moreover, the cell indices must be unique within the subset.
%           Repeated indices are not supported.
%
% RETURNS:
%   subG  - Resulting subgrid.  The cells of 'subG' are ordered according
%           to increasing numeric index in 'cells'.
%
%   gc, gf, gn
%         - Global cells, faces, and nodes referenced by sub-grid 'subG'.
%           Specifically::
%
%              gc(i) == global cell index (from 'G') of subG cell 'i'.
%              gf(i) == global face index (from 'G') of subG face 'i'.
%              gn(i) == global node index (from 'G') of subG node 'i'.
%
% NOTE:
%   The return value `gcells` is equal to `unique(cells(:))`.  Repeated
%   values in `cells` are silently ignored.
%
% EXAMPLE:
%   G    = cartGrid([3, 5, 7]);     [I, J, K] = ndgrid(2, 2:4, 3:5);
%   subG = extractSubgrid(G, sub2ind(G.cartDims, I(:), J(:), K(:)));
%   plotGrid(subG, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.25)
%   view(-35,20), camlight
%
% SEE ALSO:
%   `grid_structure`

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

   % We don't support multi-component grids...
   if numel(G) > 1
      error(msgid('Grid:MultiComponent'), ...
            'Cannot extract sub-grid from more than one grid at a time.');
   end

   if islogical(c)
      if numel(c) == G.cells.num
         c = find(c);
      else
         error(msgid('Mask:WrongSize'), ...
              ['Size of logical mask must match number of grid ', ...
               'cells.  Expected ''%d'', but got ''%d''.'], ...
               G.cells.num, numel(c));
      end
   end

   if isfield(G.cells, 'numFaces')
      error(msgid('Grid:Cells:numFaces'), ...
            'Grids with ''cells.numFaces'' are no longer supported.');
   end

   if isfield(G.faces, 'numNodes')
      error(msgid('Grid:Faces:numNodes'), ...
            'Grids with ''faces.numNodes'' are no longer supported.');
   end

   ix = @(p,i) mcolon(double(p(i)), double(p(i+1)) - 1);

   % Sort the subset cells, c, in ascending order.  Add one for outside.
   cells      = false([G.cells.num + 1, 1]);
   cells(c+1) = true;

   % replace c by sorted version of c
   c = find(cells)-1;

   % Extract faces connected to 'c'.
   faces = false([G.faces.num, 1]);
   faces(G.cells.faces(ix(G.cells.facePos, c))) = true;

   % Extract nodes connected to 'faces' (in cell subset).
   nodes = false([G.nodes.num, 1]);
   ff    = find(faces);
   nodes(G.faces.nodes(ix(G.faces.nodePos, ff))) = true;

   H.nodes.coords    = G.nodes.coords(nodes, :);

   numFaces = diff(G.cells.facePos);
   numFaces = numFaces(cells(2:end));
   numNodes = diff(G.faces.nodePos);
   numNodes = numNodes(faces);

   H.faces.neighbors = G.faces.neighbors(faces, :);

   pos = @(n) cumsum([1; double(reshape(n, [], 1))]);
   H.cells.facePos   = pos(numFaces);
   H.faces.nodePos   = pos(numNodes);

   H.cells.faces       = G.cells.faces(ix(G.cells.facePos, c ), :);
   H.faces.nodes       = G.faces.nodes(ix(G.faces.nodePos, ff)   );

   % Renumbering of grid entities (cells,faces,nodes) in subgrid.
   mc = zeros([G.cells.num+1, 1]);   mc(cells) = 1:numel(c);
   fc = zeros([G.faces.num  , 1]);   fc(faces) = 1:numel(ff);
   nc = zeros([G.nodes.num  , 1]);   nc(nodes) = 1:sum(nodes);

   H.faces.neighbors  = mc(H.faces.neighbors+1);
   H.cells.faces(:,1) = fc(H.cells.faces(:,1));
   H.faces.nodes      = nc(H.faces.nodes);

   H.faces.num = numel(H.faces.nodePos) - 1;
   H.cells.num = numel(H.cells.facePos) - 1;
   H.nodes.num = size(H.nodes.coords, 1);

   if isfield(G, 'cartDims')
      H.cartDims = nan(size(G.cartDims));
   end

   if isfield(G.cells, 'indexMap')
      H.cells.indexMap = G.cells.indexMap(cells(2:end));
      if isfield(G, 'cartDims')
         H.cartDims = G.cartDims;
      end
   end

   if isfield(G.cells, 'cpnodes')
      H.cells.cpnodes = reshape(nc(G.cells.cpnodes(cells(2:end), :)), ...
                                H.cells.num, []);
   end
   
   % Record history.
   assert (isfield(G, 'type'), ...
           'Every grid must record its origin in the field ''type''.');
   H.type    = [G.type, { mfilename }];

   assert (isfield(G, 'griddim'), ...
          ['Every grid must declare its intrinsic dimension ', ...
           'in the field ''griddim''.']);
   H.griddim = G.griddim;

   % Preserve 'computeGeometry' fields, if present.
   H = preserve_geometry(H, G, cells, faces);

   gcells = find(mc) - 1;
   gfaces = find(fc);
   gnodes = find(nc);

   % store local-to-global maps
   H.cells.global = gcells;
   H.faces.global = gfaces;
   H.nodes.global = gnodes;
end

%--------------------------------------------------------------------------

function H = preserve_geometry(H, G, cells, faces)
   if isfield(G.cells, 'volumes')
      H.cells.volumes = G.cells.volumes(cells(2:end));
   end

   if isfield(G.cells, 'centroids')
      H.cells.centroids = G.cells.centroids(cells(2:end), :);
   end

   if isfield(G.faces, 'areas')
      H.faces.areas = G.faces.areas(faces);
   end

   if isfield(G.faces, 'tag')
      H.faces.tag = G.faces.tag(faces);
   end

   if isfield(G.faces, 'normals')
      H.faces.normals = G.faces.normals(faces, :);
   end

   if isfield(G.faces, 'centroids')
      H.faces.centroids = G.faces.centroids(faces, :);
   end
end
