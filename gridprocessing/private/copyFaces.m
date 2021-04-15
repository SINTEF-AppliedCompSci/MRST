function [fnodes, nnodes, neigh, tags] = copyFaces(G, f)
%Copy face data for faces F from grid structure
%
% SYNOPSIS:
%   [fnodes, nnodes, neigh, tags] = copyFaces(G, f)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
%   f       - Face number to be copied.
%
% RETURNS:
%   fnodes  - face nodes
%   nnodes  - number of nodes for each face
%   neigh   - cell neighbors of each face
%   tags    - cellFaces tags for each half-face associated with face.
%
% COMMENTS:
%
%
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


   %   ind    = false(G.faces.num, 1); ind(f)=true;
   %   fnodes = G.faces.nodes(rldecode(ind(:), double(G.faces.numNodes)));

   p1 = G.faces.nodePos(f);
   p2 = G.faces.nodePos(f+1)-1;
   fnodes = G.faces.nodes(mcolon(p1,p2));
   neigh  = G.faces.neighbors(f,:);
   nnodes = diff(G.faces.nodePos);
   nnodes = nnodes(f);

   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
   col    = 1+(cellNo~=G.faces.neighbors(G.cells.faces(:,1), 1));
   if size(G.cells.faces,2)>1,
      tags = accumarray([G.cells.faces(:,1), col], G.cells.faces(:,2));
      tags = tags(f,:);
   else
      tags =[];
   end
end
