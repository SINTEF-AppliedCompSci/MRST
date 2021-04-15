function G = processUnstructuredEGRID(egrid)
%Make grid structure from unstructured Eclipse EGRID fields.
%
% SYNOPSIS:
%   G = processUnstructuredEGRID(egrid)
%
% PARAMETERS:
%   egrid   - struct with fieldnames and data read from Eclipse EGRID data
%             structure.
% RETURNS:
%   G       - Grid structure as described by grid_structure.
%
% COMMENTS:
%
% SEE ALSO:
%   `readEclipseOutputFileUnFmt`

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


   if strcmp(egrid.GRIDUNIT.values(1), 'FEET')
      egrid.NDCOORD.values = convertFrom(egrid.NDCOORD.values, ft);
   end
   G.nodes  = struct('num',     numel(egrid.NDCOORD.values)/3, ...
                     'coords',  reshape(egrid.NDCOORD.values, 3, [])');

   i = egrid.NCELLFAC.values>0;
   G.cells  = struct('num',     numel(egrid.NCELLFAC.values(i)), ...
                     'faces',   egrid.CELLFACS.values, ...
                     'facePos', cumsum([1;egrid.NCELLFAC.values(i)]));

   G.faces  = struct('num',     numel(egrid.NFACENOD.values), ...
                     'nodes',   egrid.FACENODS.values, ...
                     'nodePos', cumsum([1;egrid.NFACENOD.values]));

   G.faces.neighbors = neighborsFromCellFaces(G);

   G.faces.tag = inf(G.faces.num, 3);
   if isfield(egrid, 'FACEHING'),
      G.faces.tag(egrid.FACEHING.values(1:2:end), :) = ...
         G.nodes.coords(egrid.FACEHING.values(2:2:end), :);
   end

   G.type    = { mfilename };
   G.griddim = 3;
end

function N = neighborsFromCellFaces(G)
%Make num-faces x 2 cell neighbors array from other grid fields.
%
% SYNOPSIS:
%   G = neighborsFromCellFaces(G)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
% RETURNS:
%   N       - An n x 2 array of cell numbers.
%
% COMMENTS:
%
% SEE ALSO:
%

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

   N(cellfaces(isboundary), 1) = cellNo(isboundary);
   N = findNormalDirections(G, N);
end

function N = findNormalDirections(G, N)
% Fix neighbour array to conform to sign convention for faces.
%
% SYNOPSIS:
%   N = findNormalDirections(G, N)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
%   N       - n x 2 cell neighbors array.
%
% RETURNS:
%   N       - n x 2 cell neighbors array.

   assert(size(G.nodes.coords, 2)==3);

   % Assume convex faces.   Compute average of node coordinates.
   faceno = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2) .';
%%{
   fcenters = zeros(G.faces.num, size(G.nodes.coords, 2));
   for d = 1:size(G.nodes.coords, 2),
      fcenters(:,d) = accumarray(faceno, G.nodes.coords(G.faces.nodes, d))./...
         accumarray(faceno,1);
   end
%}
%{
   fcenters = sparse(faceno, G.faces.nodes(:,1), 1) ...
              * [G.nodes.coords, ones([G.nodes.num, 1])];
   fcenters = bsxfun(@rdivide, fcenters(:, 1:end-1), fcenters(:,end));
%}

   % Assume convex cells.
   cellno = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
%%{
   ccenters = zeros(G.cells.num, size(G.nodes.coords, 2));
   for d = 1:size(G.nodes.coords, 2),
      ccenters(:,d) = accumarray(cellno, fcenters(G.cells.faces(:,1), d))./...
         accumarray(cellno, 1);
   end
%}
%{
   ccenters = sparse(cellno, G.cells.faces(:,1), 1) ...
              * [fcenters, ones([G.faces.num, 1])];
   ccenters = bsxfun(@rdivide, ccenters(:, 1:end-1), ccenters(:,end));
%}

   % Compute triple product v1 x v2 · v3 of vectors v1 = fc-cc, v2 = n1-fc,
   % and v3 = n2-n1 --- cc and fc being cell centers and face centers, n1
   % and n2 being the first and second node of the face.  Triple product
   % should be positive for half-faces with positive sign.

   assert (all(diff(G.faces.nodePos) >= 3), ...
           'All faces must have at least three nodes in 3D.');

   n1 = G.nodes.coords(G.faces.nodes(G.faces.nodePos(1:end-1)), :);
   n2 = G.nodes.coords(G.faces.nodes(G.faces.nodePos(1:end-1)+1), :);
   v1 = fcenters(G.cells.faces(:,1), :) - ccenters(cellno,:);
   v2 = n1(G.cells.faces(:,1), :) - fcenters(G.cells.faces(:,1), :);
   v3 = n2(G.cells.faces(:,1),:) - n1(G.cells.faces(:,1),:);

   a   = sum(cross(v1,v2).*v3, 2);
   sgn = -1.0+2.0*(N(G.cells.faces(:,1), 1) == cellno);
   i   = accumarray(G.cells.faces(:,1), a.*sgn)<0;
   N(i,:)= N(i,[2,1]);
end
