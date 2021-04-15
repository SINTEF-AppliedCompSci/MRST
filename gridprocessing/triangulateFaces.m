function G = triangulateFaces(G, f)
%Split face f in grid G into subfaces.
%
% DESCRIPTION:
%   For each face listed in F, face center as compute average of face node
%   coordinates.  Use face center to triangulate face.  remove original
%   face and add triangular faces.
%
%   If no face list F is given, all faces in G are triangulated.  This
%   should make a grid with curved faces polyhedral, i.e., a grid with onlu
%   plane faces.
%
% SYNOPSIS:
%   G = triangulateFaces(G)
%   G = triangulateFaces(G, f)
%
% PARAMETERS:
%   G        - Grid
%
%   f        - indices to faces that are to be split. Default is all faces.
%
% RETURNS:
%   G        - Grid.
%
% SEE ALSO:
%   `grid_structure`, `computeGeometry`

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


   if nargin == 1,
      f = (1:G.faces.num)';
   end

   [G, newnodes] = addMidpointNodes(G, f);

   % construct new faces in terms of (newfacenodes, newneighbors)
   ix  = mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1);
   nix = G.faces.nodes(ix);
   j   = double(G.faces.nodePos(f+1)-G.faces.nodePos(f));
   d = size(G.nodes.coords,2);
   if d == 2,
      t = [nix(1:2:end), newnodes,newnodes, nix(2:2:end)];
   else
      pos = cumsum([0;j]);
      t   = [nix, [nix(2:end);0], rldecode(newnodes, j)];
      t(pos(2:end),2) = t(1+pos(1:end-1),1);
   end

   newfacenodes = reshape(t',[], 1);
   newneighbors = rldecode(G.faces.neighbors(f,:), j);
   newnumnodes  = repmat(d, [size(newneighbors,1), 1]);

   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
   col    = 1+(cellNo~=G.faces.neighbors(G.cells.faces(:,1), 1));
   hftags = accumarray([G.cells.faces(:,1), col], G.cells.faces(:,2));
   hftags = rldecode(hftags(f,:), j);

   G      = removeFaces(G, f);
   G      = addFaces(G, newfacenodes, newnumnodes, newneighbors, hftags);

   G.type = [G.type, { mfilename }];
end

function [H, newnodes] = addMidpointNodes(G, f)
   d = size(G.nodes.coords, 2);

   ix     = mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1);
   faceNo = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2) .';
   nodes  = zeros(numel(f), d);
   for i=1:d
      v = accumarray(faceNo(ix), G.nodes.coords(G.faces.nodes(ix), i));
      nodes(:,i) = v(f);
   end

   denominator    = accumarray(faceNo(ix), 1);
   nodes          = bsxfun(@rdivide, nodes, denominator(f));

   newnodes       = double(G.nodes.num) +(1:size(nodes,1))';
   G.nodes.coords = [G.nodes.coords; nodes];
   G.nodes.num    = size(G.nodes.coords,1);
   H = G;
end

