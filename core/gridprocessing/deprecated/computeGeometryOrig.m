function G = computeGeometryOrig(G, varargin)
%Compute geometry of grid.
%
% SYNOPSIS:
%   G = computeGeometryOrig(G)
%   G = computeGeometryOrig(G, 'pn', pv, ...)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
%   'pn'/pv - List of 'key'/value pairs for supplying optional parameters.
%             The supported options are
%               - Verbose -- Whether or not to display informational
%                            messages throughout the computational process.
%                            Logical.  Default value: Verbose = false
%                            (don't display any informational messages).
% OPTIONAL PARAMETERS:
%
%   'hingenodes'  - A struct with fields 'faces' and 'nodes'.  A hinge node
%                   is an extra center node for a face, that is used to
%                   triangulate the face geometry.  For each face number F
%                   in 'faces' there is a row in 'nodes' which holds the
%                   node coordinate for the hinge node belonging to face F.
%
% RETURNS:
%   G - Grid structure with added fields:
%         - cells
%             - volumes   -- A G.cells.num-by-1 array of cell volumes.
%
%             - centroids -- A G.cells.num-by-SIZE(G.nodes.coords, 2) array
%                            of (approximate) cell centroids.
%
%         - faces
%             - areas     -- A G.faces.num-by-1 array of face areas.
%
%             - normals   -- A G.faces.num-by-G.griddim array
%                            of face normals.
%
%             - centroids -- A G.faces.num-by-SIZE(G.nodes.coords, 2) array
%                            of (approximate) face centroids.
%
% COMMENTS:
%   Individual face normals have length (i.e., Euclidian norm) equal to
%   the corresponding face areas.  In other words, subject to numerical
%   round-off, the identity
%
%         NORM(G.faces.normals(i,:), 2) == G.faces.areas(i)
%
%   holds for all faces i=1:G.faces.num .
%
%   In three space dimensions, i.e., when G.griddim == 3,
%   function 'computeGeometry' assumes that the nodes on a given face, f,
%   are ordered such that the face normal on f is directed from cell
%   G.faces.neighbors(f,1) to cell G.faces.neighbors(f,2).
%
% SEE ALSO:
%   `grid_structure`.

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


%% Setup
assert(size(G.faces.nodes, 2)==1);
opt     = struct('verbose',              mrstVerbose, ...
                 'findNeighbors',        false,       ...
                 'hingenodes',           []);
opt     = merge_options(opt, varargin{:});

assert(isempty(opt.hingenodes) || G.griddim == 3, ...
   'Hinge nodes are only supported for 3D grids.');

numC    = G.cells.num;
numF    = G.faces.num;

%% Possibly find neighbors
if opt.findNeighbors,
   G.faces.neighbors = findNeighbors(G);
   G = findNormalDirections(G);
else
   if ~isfield(G.faces, 'neighbors'),
      warning(msgid('GridType:Incomplete'), ...
         ['No field ''faces.neighbors'' found. ',...
         'Adding plausible values... proceed with caution!']);
      G.faces.neighbors = findNeighbors(G);
      G = findNormalDirections(G);
   end
end

%% Main part
if G.griddim == 3,
   assert(size(G.nodes.coords,2) == 3);
   faceNo  = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2) .';
   p       = G.faces.nodePos;
   next    = (2:size(G.faces.nodes, 1)+1) .';
   next(p(2 : end) - 1) = p(1 : end-1);

   %%
   % Divide each face into sub-triangles all having one node as
   % pCenter = sum(nodes) /#nodes. Compute area-weighted normals,
   % and add to obtain approx face-normals. Compute resulting areas
   % and centroids.
   dispif(opt.verbose, 'Computing normals, areas, and centroids...\t');
   t0 = ticif (opt.verbose);

   llE = length(G.faces.nodes);
   localEdge2Face = sparse(1 : llE, faceNo, 1, llE, numF);

   pCenters     = bsxfun(@rdivide, ...
                         localEdge2Face.' * G.nodes.coords(G.faces.nodes,:), ...
                         diff(double(G.faces.nodePos)));
   pCenters     = localEdge2Face * pCenters;

   % Use hinge nodes for selected faces if present.
   if ~isempty(opt.hingenodes),
      ix = mcolon(G.faces.nodePos(opt.hingenodes.faces), ...
                  G.faces.nodePos(opt.hingenodes.faces+1)-1);
      pCenters(ix, :) = rldecode(opt.hingenodes.nodes, ...
                        G.faces.nodePos(opt.hingenodes.faces+1) - ...
                        G.faces.nodePos(opt.hingenodes.faces));
   end

   subNormals   = cross(G.nodes.coords(G.faces.nodes(next),:) - ...
                        G.nodes.coords(G.faces.nodes,:), ...
                        pCenters - G.nodes.coords(G.faces.nodes,:)) ./ 2;
   subAreas     = sqrt(sum(subNormals .^ 2, 2));
   subCentroids = (G.nodes.coords(G.faces.nodes,:) + ...
                   G.nodes.coords(G.faces.nodes(next),:) + pCenters) ./ 3;
   clear llE faceNo

   faceNormals    = localEdge2Face.' * subNormals;
   faceAreas      = localEdge2Face.' * subAreas;
   subNormalSigns = sign(sum(subNormals .* (localEdge2Face * faceNormals), 2));
   faceCentroids  = bsxfun(@rdivide,                                 ...
                           localEdge2Face.' * ...
                           bsxfun(@times, subAreas, subCentroids), ...
                           faceAreas);

   % Computation above does not make sense for faces with zero area
   i=find(faceAreas==0);
   if numel(i)>0,
      warning(msgid('computeGeometry:faceAreas'), ...
         [' Faces with zero area detected. Such faces should be ', ...
          'removed before calling ', mfilename]);
       faceCentroids(i,:) = pCenters(i,:);
   end

   clear subAreas pCenters

   tocif(opt.verbose, t0)

   %%
   % Divide each cell into sub-tetrahedra according to sub-triangles above,
   % all having one node as cCenter = sum(faceCentroids) / #faceCentroids.

   dispif(opt.verbose, 'Computing cell volumes and centroids...\t\t');
   t0 = ticif (opt.verbose);

   cellVolumes   = zeros([numC, 1]);
   cellCentroids = zeros([numC, 3]);

   lastInx = 0;
   for c = 1 : numC,
      nF  = double(G.cells.facePos(c+1)- G.cells.facePos(c));
      inx = (1 : nF) + lastInx;

      faces        = G.cells.faces(inx,1);
      [triE, triF] = find(localEdge2Face(:,faces));

      fCentroids = faceCentroids(faces,:);
      cCenter    = sum(fCentroids) ./ double(nF);

      relSubC    = bsxfun(@minus, subCentroids(triE,:), cCenter);

      % The normal of a face f is directed from cell G.faces.neighbors(f,1)
      % to cell G.faces.neighbors(f,2).   If cell c is in the second column
      % for face f, then the nomal must be multiplied by -1 to be an outer
      % normal.
      orientation = 2*double(G.faces.neighbors(G.cells.faces(inx,1), 1) == c)-1;

      outNormals = bsxfun(@times,             ...
                          subNormals(triE,:), ...
                          subNormalSigns(triE) .* orientation(triF));

      tVolumes   = (1/3) * sum(relSubC .* outNormals, 2);
      tCentroids = (3/4) * relSubC;

      volume      = sum(tVolumes);
      relCentroid = (tVolumes' * tCentroids) ./ volume;
      centroid    = relCentroid + cCenter;

      cellVolumes(c)     = volume;
      cellCentroids(c,:) = centroid;

      lastInx = lastInx + nF;
   end
   tocif(opt.verbose, t0)

elseif (G.griddim ==2) && (size(G.nodes.coords,2)==2)

  dispif(opt.verbose, 'Computing normals, areas, and centroids...\t');
  t0 = ticif (opt.verbose);

  edges     = reshape(G.faces.nodes,2,[])';
  edgeLength    = G.nodes.coords(edges(:,2),:)-G.nodes.coords(edges(:,1),:);
  faceAreas     = sqrt(sum(edgeLength.*edgeLength,2));
  faceCentroids = (G.nodes.coords(edges(:,2),:)+ ...
                   G.nodes.coords(edges(:,1),:))/2;
  faceNormals   = [edgeLength(:,2),-edgeLength(:,1)];

  tocif(opt.verbose, t0)
  dispif(opt.verbose, 'Computing cell volumes and centroids...\t\t');
  t0 = ticif (opt.verbose);


  numfaces=diff(G.cells.facePos);

  cellno        = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
  cellEdges     = edges(G.cells.faces(:,1),:);
  r             = G.faces.neighbors(G.cells.faces(:,1), 2) == cellno;
  cellEdges(r,:)= cellEdges(r, [2,1]);

  cCenter       = zeros(G.cells.num, 2);
  cCenter(:,1)  = accumarray(cellno, faceCentroids(G.cells.faces(:,1), 1));
  cCenter(:,2)  = accumarray(cellno, faceCentroids(G.cells.faces(:,1), 2));
  cCenter       = bsxfun(@rdivide, cCenter, double(numfaces));

  a             = G.nodes.coords(cellEdges(:,1),:) - cCenter(cellno,:);
  b             = G.nodes.coords(cellEdges(:,2),:) - cCenter(cellno,:);
  subArea       = 0.5*(a(:,1).*b(:,2)-a(:,2).*b(:,1));

  subCentroid   = bsxfun(@plus, cCenter(cellno,:),  2*faceCentroids(G.cells.faces(:,1),:))/3;
  cellVolumes   = accumarray(cellno, subArea);

  cellCentroids = zeros(G.cells.num, 2);
  cellCentroids(:,1) = accumarray(cellno, bsxfun(@times, subArea, subCentroid(:,1)));
  cellCentroids(:,2) = accumarray(cellno, bsxfun(@times, subArea, subCentroid(:,2)));
  cellCentroids = bsxfun(@rdivide, cellCentroids, cellVolumes);

  tocif(opt.verbose, t0)

elseif (G.griddim ==2) && (size(G.nodes.coords,2)==3)
  dispif(opt.verbose, ...
     'Experimental implementation only available for surface grids\n');

  dispif(opt.verbose, 'Computing normals, areas, and centroids...\t');
  t0 = ticif (opt.verbose);

  edges     = reshape(G.faces.nodes,2,[])';
  edgeLength    = G.nodes.coords(edges(:,2),:)-G.nodes.coords(edges(:,1),:);
  faceAreas     = sqrt(sum(edgeLength.*edgeLength,2));
  faceCentroids = (G.nodes.coords(edges(:,2),:)+ ...
                   G.nodes.coords(edges(:,1),:))/2;
  faceNormals   = [edgeLength(:,2),-edgeLength(:,1)];

  tocif(opt.verbose, t0)
  dispif(opt.verbose, 'Computing cell volumes and centroids...\t\t');
  t0 = ticif (opt.verbose);

  numfaces=diff(G.cells.facePos);

  cellno        = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
  cellEdges     = edges(G.cells.faces(:,1),:);
  r             = G.faces.neighbors(G.cells.faces(:,1), 2) == cellno;
  cellEdges(r,:)= cellEdges(r, [2,1]);

  cCenter       = zeros(G.cells.num, 3);
  cCenter(:,1)  = accumarray(cellno, faceCentroids(G.cells.faces(:,1), 1));
  cCenter(:,2)  = accumarray(cellno, faceCentroids(G.cells.faces(:,1), 2));
  cCenter(:,3)  = accumarray(cellno, faceCentroids(G.cells.faces(:,1), 3));
  cCenter       = bsxfun(@rdivide, cCenter, double(numfaces));

  a             = G.nodes.coords(cellEdges(:,1),:) - cCenter(cellno,:);
  b             = G.nodes.coords(cellEdges(:,2),:) - cCenter(cellno,:);
  subArea       = sqrt(sum(cross(a,b,2).^2,2))*0.5;

  subCentroid   = bsxfun(@plus, cCenter(cellno,:), ...
                         2*faceCentroids(G.cells.faces(:,1),:))/3;
  cellVolumes   = accumarray(cellno, subArea);

  cellCentroids = zeros(G.cells.num, 3);
  cellCentroids(:,1) = accumarray(cellno, bsxfun(@times, subArea, subCentroid(:,1)));
  cellCentroids(:,2) = accumarray(cellno, bsxfun(@times, subArea, subCentroid(:,2)));
  cellCentroids(:,3) = accumarray(cellno, bsxfun(@times, subArea, subCentroid(:,3)));
  cellCentroids = bsxfun(@rdivide, cellCentroids, cellVolumes);

  tocif(opt.verbose, t0)
else
  assert(false);
end

%% Update grid
G.faces.areas     = faceAreas;
G.faces.normals   = faceNormals;
G.faces.centroids = faceCentroids;

G.cells.volumes   = cellVolumes;
G.cells.centroids = cellCentroids;

if ~isfield(G, 'type'),
   warning(msgid('GridType:Unknown'),                            ...
          ['Input grid has no known type. ',                     ...
           'I''ll assume it arose from the primordial soup...']);
   G.type = { 'Primordial Soup' };
end

G.type = [G.type, { mfilename }];
end

%--------------------------------------------------------------------------

function N = findNeighbors(G)

   % Internal faces
   cellNo         = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   [cellfaces, j] = sort(G.cells.faces(:,1));
   cellNo         = cellNo(j);
   hf             = find(cellfaces(1:end-1) == cellfaces(2:end));

   N                       = zeros(G.faces.num, 1);
   N(cellfaces(hf,1),   1) = cellNo(hf);
   N(cellfaces(hf + 1), 2) = cellNo(hf + 1);

   % Boundary faces
   isboundary         = true(numel(cellNo), 1);
   isboundary(hf)     = false;
   isboundary(hf + 1) = false;
   %hf                 = find(isboundary);
   N(cellfaces(isboundary), 1) = cellNo(isboundary);
end

%--------------------------------------------------------------------------

function G = findNormalDirections(G)
   assert (all([G.griddim, size(G.nodes.coords, 2)] == 3), ...
          ['Detecting neighbourship based on normal directions ', ...
           'is only supported in 3D grids.']);

   % Assume convex faces.   Compute average of node coordinates.
   fcenters = ...
      averageCoordinates(diff(G.faces.nodePos), ...
                         G.nodes.coords(G.faces.nodes, :));

   % Assume convex cells.   Compute average of face centre coordinates.
   [ccenters, cellno] = ...
      averageCoordinates(diff(G.cells.facePos), ...
                         fcenters(G.cells.faces(:,1), :));

   % Compute triple product v1 x v2 · v3 of vectors v1 = fc-cc, v2 = n1-fc,
   % and v3 = n2-n1 --- cc and fc being cell centres and face centres, n1
   % and n2 being the first and second node of the face.  Triple product
   % should be positive for half-faces with positive sign.

   n1 = G.nodes.coords(G.faces.nodes(G.faces.nodePos(1:end-1)    ), :);
   n2 = G.nodes.coords(G.faces.nodes(G.faces.nodePos(1:end-1) + 1), :);

   v1 = fcenters(G.cells.faces(:,1), :) - ccenters(cellno, :);
   v2 = n1(G.cells.faces(:,1), :) - fcenters(G.cells.faces(:,1), :);
   v3 = n2(G.cells.faces(:,1), :) - n1(G.cells.faces(:,1), :);

   a   = sum(cross(v1, v2) .* v3, 2);
   sgn = 2*(G.faces.neighbors(G.cells.faces(:,1), 1) == cellno) - 1;

   i   = accumarray(G.cells.faces(:,1), a .* sgn) < 0;
   G.faces.neighbors(i, :) = G.faces.neighbors(i, [2, 1]);
end

%--------------------------------------------------------------------------

function [c, no] = averageCoordinates(n, c)
   no = rldecode(1 : numel(n), n, 2) .';
   c  = sparse(no, 1 : numel(no), 1) * [ c, ones([size(c, 1), 1]) ];
   c  = bsxfun(@rdivide, c(:, 1 : end - 1), c(:, end));
end
