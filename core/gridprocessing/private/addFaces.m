function [H, newfaces] = addFaces (G, facenodes, numnodes, neighbors, tags)
%Add faces F from grid structure
%
% SYNOPSIS:
%   H             = addFaces(G, facenodes, numnodes, neighbors, tags)
%   [H, newfaces] = addFaces(...)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
%   fnodes  - face nodes
%
%   nnodes  - number of nodes for each face
%   neigh   - cell neighbors of each face
%   tags    - cellFaces tags for each half-face associated with face.
%
% RETURNS:
%   G       - Grid structure where the following fields have been modified:
%
%                  G.faces.nodes
%
%                  G.cells.faces
%                  G.cells.facePos
%                  G.cells.numFaces
%
%                  G.faces.neighbors
%                  G.faces.numNodes
%                  G.faces.nodePos
%                  G.faces.tag
%                  G.faces.num
%
%  newfaces - Numbering of new faces.
%
% COMMENTS:
%   Cells may not be closed polygons/polyhedra after this call.

% BUG
%
% o Does not preserve orientation in 2D

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

   assert(numel(numnodes) == size(neighbors, 1));
   assert(sum(numnodes) == numel(facenodes));

   G.faces.neighbors = [G.faces.neighbors; neighbors];
   %G.faces.tag       = [G.faces.tag; tags];
   G.faces.nodes     = [G.faces.nodes; facenodes];

   numf        = numel(numnodes);
   newfaces    = G.faces.num + (1 : numf)';
   G.faces.num = double(G.faces.num + numf);
   clear numf

   if isfield(G.faces, 'numNodes')
      G.faces.numNodes = [G.faces.numNodes; numnodes];
      G.faces.nodePos = cumsum([1; double(G.faces.numNodes)]);
   else
      G.faces.nodePos = cumsum([1; double(diff(G.faces.nodePos)); ...
                                double(numnodes)]);
   end

   faces = repmat(newfaces, [2,1]);
   k     = neighbors(:) ~= 0;

   if nargin == 5
      tags = tags(k);
   end

   faces     = faces(k);
   neighbors = neighbors(k);

   if nargin == 4
      [G.cells.faces,G.cells.facePos] = ...
         insertInPackedData(G.cells.facePos, G.cells.faces, ...
                            neighbors(:), faces(:));
   else
      [G.cells.faces,G.cells.facePos] = ...
         insertInPackedData(G.cells.facePos, G.cells.faces, ...
                            neighbors(:), [faces(:), tags(:)]);
   end

   % keep for the time being
   if isfield(G.faces, 'tag')
      G.faces.tag = [G.faces.tag; zeros(numel(numnodes), size(G.faces.tag, 2))];
   end

   H = G;
end
