function [H, facemap] = removeFaces(G, f)
%Remove faces F from grid structure
%
% SYNOPSIS:
%   H = removeFaces(G, f)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
%   f       - Face number to be removed.
%
% RETURNS:
%   G       - Grid structure where the followinf fields have been modified:
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
%   facemap - Index map that allow conversion from old face numbers to new
%             face numbers.
%
% COMMENTS:
%  Cells may not be closed polygons/polyhedra after this call.

% SEE ALSO:
%

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


   ind     = false(G.faces.num, 1); ind(f)=true;
   facemap = mapExcluding(ind);

   % remove nodes associated with faces f in G.faces.nodes
   numNodes = diff(G.faces.nodePos);
   G.faces.nodes(rldecode(ind(:), double(numNodes))) = [];

   % remove and renumber faces in cellFaces
   [G.cells.faces, G.cells.facePos] = ...
      removeFromPackedData(G.cells.facePos, G.cells.faces, f);

   G.cells.faces(:,1) = facemap(G.cells.faces(:,1));

   G.faces.neighbors(ind,:) = [];
   numNodes (ind,:) = [];
   G.faces.nodePos = cumsum([1; double(numNodes)]);
   if isfield(G.faces, 'numNodes'),
      G.faces = rmfield(G.faces, 'numNodes');
   end
   if isfield(G.faces, 'tag'),
      G.faces.tag (ind,:) = [];
   end

   G.faces.num = G.faces.num - sum(ind);

   H = G;
end

%--------------------------------------------------------------------------

function m = mapExcluding(indices)
   n            = numel(indices);
   ind          = ones(n,1);
   ind(indices) = 0;
   m            = cumsum(ind);
   m(indices)   = 0;
end
