function G = makePlanarGrid(H, faces, sgn)
%Construct 2D surface grid from faces of 3D grid.
%
% SYNOPSIS:
%   g = makePlanarGrid(G, faces, sgn)
%
% PARAMETERS:
%   G      - Valid grid definition containing connectivity, cell
%            geometry, face geometry and unique nodes.
%
%   faces  - List of faces in original 3D grid G that are cells in g.
%
%   sgn    - Sign of face in 3D grid associated with each cell in g.
%
% RETURNS:
%   g      - Valid 2D grid.
%
% EXAMPLE:
%   H = cartGrid([2,3,4]);
%   G = makePlanarGrid(H, (1:3:36)', ones(12, 1));
%   plotGrid(G); plotGrid(H, 'facecolor','none');view(3);
%
%
% NOTE:
%   All node coordinates in `H.nodes` are copied to `G.nodes`, and
%   `G.faces.nodes` refer to this numbering of nodes.
%
% SEE ALSO:
%   `makePlanarGrid`

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

% Written by Jostein Natvig, SINTEF Applied Mathematics.

% Make 2D grid embedded in 3D from 3D faces
   if nargin < 3, sgn = ones(numel(faces), 1); end
   pos = [H.faces.nodePos(faces), H.faces.nodePos(faces+1)-1];
   pos(sgn<0,:) = pos(sgn<0,[2,1]);
   fnodes     = double(H.faces.nodes(mcolon(pos(:,1), pos(:,2), sgn)));
   nnodes     = double(abs(diff(pos, 1, 2))+1);
   pos        = cumsum([1; nnodes]);
   cellno     = rldecode((1:numel(nnodes))', nnodes);
   edges      = [fnodes, fnodes(rot(pos, 1))];

   % Identify unique edges
   [e,i]      = sort(edges, 2);
   [j,j]      = sortrows(e);
   k(j)       = 1:numel(j);
   [edges, n] = rlencode(e(j,:));
   edgenum    = rldecode((1:size(edges,1))', n);
   edgenum    = edgenum(k);
   edgesign   = i(:,1);

   % Identify neigbors assuming two per edge...
   p          = sub2ind(size(edges), edgenum, edgesign);
   neigh      = zeros(size(edges));
   neigh(p)   = cellno;

   G.faces.num        = size(edges, 1);
   G.faces.nodes      = reshape(edges', [], 1);
   G.faces.nodePos    = (1:2:2*G.faces.num+1)';
   G.faces.neighbors  = neigh;
   G.cells.num        = numel(nnodes);
   G.cells.facePos    = pos;
   G.cells.faces(:,1) = edgenum;
   G.cells.faces(:,2) = 0;
   G.nodes            = H.nodes;
   G.dim              = 2;
   G.griddim          = 2;
end

% Useful operation on position vectors: Rotation modulo section size
function ix= rot(pos, offset)
   num    = diff(pos);
   offset = mod(offset, num); % net offset
   ix     = zeros(max(pos)-1, 1);

   ix(mcolon(pos(1:end-1), pos(1:end-1)+num-offset-1)) = ...
      mcolon(pos(1:end-1)+offset, pos(2:end)-1);
   ix(mcolon(pos(1:end-1)+num-offset, pos(2:end)-1)) = ...
      mcolon(pos(1:end-1), pos(2:end)-1-num+offset);
end
