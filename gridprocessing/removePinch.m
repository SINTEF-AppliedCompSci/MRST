function G = removePinch(G, tol)
%Uniquify nodes, remove pinched faces and cells.
%
% SYNOPSIS:
%   H = removePinch(G)
%   H = removePinch(G, tol)
%
% PARAMETERS:
%   G       - Grid structure as described by `grid_structure`.
%   tol     - Absolute tolerance to distinguish neighbouring points
%
% RETURNS:
%   G       - Grid structure where duplicate nodes have been removed in
%             `G.nodes.coords` and duplicate node numbers are removed from
%             `G.faces.nodes`.  Faces with fewer than 3 nodes and cells with
%             fewer than `G.griddim+1` faces are also removed to avoid zero
%             areas and zero volumes subsequently.
%
% SEE ALSO:
%  `removeCells`, `extractSubgrid`

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

% Written by Jostein R. Natvig, SINTEF Applied Mathematics.

   % Uniquify nodes
   [G.nodes.coords, i, j] = unique(G.nodes.coords, 'rows');

   % Map face nodes
   G.faces.nodes    = j(G.faces.nodes);

   % Remove nodes with small difference
   if nargin == 2
      d = [inf; sqrt(sum(diff(G.nodes.coords,1) .^ 2, 2))];
      I = d < tol;
      G.nodes.coords = G.nodes.coords(~I,:);
      J = ones(size(I));
      J(I) = 0;  J = cumsum(J);
      G.faces.nodes = J(G.faces.nodes);
   end
   G.nodes.num = size(G.nodes.coords, 1);

   % remove repeated node numbers in faces (stored in adjacent positions
   faceno = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2)';
   tmp    = rlencode([faceno, G.faces.nodes]);
   [n,n]  = rlencode(tmp(:,1));

   % remove nodes that coincide without being stored in adjacent positions
   pos    = cumsum([1;n]);
   ix     = tmp(pos(1:end-1), 2) == tmp(pos(2:end)-1, 2) & ...
               (pos(1:end-1)     ~=     pos(2:end)-1);
   tmp(pos(ix),:)  = [];

   G.faces.nodes   = tmp(:, 2);
   [fn,n]          = rlencode(tmp(:,1));
   G.faces.nodePos = cumsum([1; n]);
   G.faces.num     = numel(G.faces.nodePos)-1;


   % remove pinched faces
   G = removeFaces(G, find(diff(G.faces.nodePos)<G.griddim));



   % Permute face nodes to canonical order
   % 1) Rotate smallest face node to first position
   p = G.faces.nodePos(1:end-1);
   N = max(diff(G.faces.nodePos));
   m = inf(G.faces.num, 1);
   v = zeros(G.faces.num, 1);
   for i=1:N,
      j    = G.faces.nodes(p) < m;
      v(j) = p(j);
      m(j) = G.faces.nodes(p(j));

      p = min([p+1, G.faces.nodePos(2:end)-1], [], 2);
   end
   offset = v-G.faces.nodePos(1:end-1);
   ix     = rot(G.faces.nodePos, offset);
   G.faces.nodes = G.faces.nodes(ix);

   % 2) flip face nodes to place 2. largest face node in second position.
   i = G.faces.nodes(G.faces.nodePos(1:end-1)+1) < ...
       G.faces.nodes(G.faces.nodePos(2:end  )-1);
   [ix1, ix2] = flip(G.faces.nodePos, i);
   G.faces.nodes(ix2) = G.faces.nodes(ix1);
   G.faces.neighbors(i,:) = G.faces.neighbors(i,[2,1]);

   % Remove cells with fewer than griddim+1 faces
   G = removeCells(G, find(diff(G.cells.facePos) < G.griddim + 1));

   % Optional: reconnect cell neighbors
   G = uniqueFaces(G);
end

function G = uniqueFaces(G)
%Find pairs of faces with different face ID but identical nodes.
%
%  Use removeInternalBoundary to merge each pair of faces to produce
%  connectivity across (totally) pinched layers.  Connectivity in partially
%  pinched zones is always retained, somehow.
%
% SYNOPSIS:
%   G = uniqueFaces(G)
%
% PARAMETERS:
%   G - Grid structure as described by grid_structure.
%
% RETURNS:
%
%   G - Modified grid structure with only unique face IDs.
%
%

   f     = find(any(G.faces.neighbors==0, 2));
   ix    = mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1);
   nodes = G.faces.nodes(ix);
   n     = double(G.faces.nodePos(f+1)-G.faces.nodePos(f));
   N     = max(n);

   % Make inf-padded array of node IDs for each "outer" face in 'f'
   v     = inf(N, numel(n));
   k     = rldecode(1:N:N*numel(n), n, 2)' + mcolon(0, n-1)';
   v(k)  = nodes;
   v     = v';

   % Find identical faces (by node ID)
   [v,i] = sortrows(v);
   [n,n] = rlencode(v);
   N     = rldecode(n,n);

   assert(all(N>0 & N<3), 'Internal error in removePinch::uniqueFaces!');

   % Remove faces, preserve faces.tags.
   F      = reshape(f(i(N==2)), 2, [])';
   if(isfield(G.faces,'tag'))
      tags   = G.faces.tag(F(:,1), :);
   else
      tags=[];
   end
   [G, f] = removeInternalBoundary(G, F);
   if(~isempty(tags))
      G.faces.tag(f,:) = tags;
   end
end


function [ix1, ix2] = flip(pos, i)
%Make index set to flip data sections in packed array
%
% SYNOPSIS:
%   H = flip(G, f)
%
% PARAMETERS:
%   pos - Indirection map to sections of data.
%
%   i   - Sections to flip.
%
% RETURNS:
%   ix1 - index set for right-hand side.
%
%   ix2 - index set for left-hand side.
%
% EXAMPLE
%
%   % Rotate nodes in first 2 faces.
%   [ix1, ix2] = flip(G.faces.nodePos, [1,2]);
%   G.faces.nodes(ix2) = G.faces.nodes(ix1);
%
   if islogical(i), i=find(i); end

   ix1 = mcolon(pos(i+1)-1, pos(i), -1);
   ix2 = mcolon(pos(i), pos(i+1)-1);
end



function ix = rot(pos, offset)
%Make index set to rotate data sections in packed array
%
% SYNOPSIS:
%   ix = rot(pos, offset)
%
% PARAMETERS:
%   pos    - Indirection map to sections of data.
%
%   offset - number of positions to rotate each section of data.
%
%   i      - which sections to rotate. Reqired: numel(i) == numel(offset).
% RETURNS:
%   ix     - index set into data settions.
%
% EXAMPLE
%
%   % Rotate nodes in all faces by 2 positions.
%   ix = rot(G.faces.nodePos, 2*ones(G.faces.num, 1));
%   G.faces.nodes = G.faces.nodes(ix);
%
   i = (1:numel(pos)-1)';
   assert(numel(i) == numel(offset));

   num    = pos(i+1)-pos(i);
   offset = mod(offset(:), num(:)); % net offset
   ix     = zeros(sum(num), 1);

   ix(mcolon(pos(i),        pos(i)+num-offset-1)) = ...
      mcolon(pos(i)+offset, pos(i+1)-1);
   ix(mcolon(pos(i)+num-offset, pos(i+1)-1)) = ...
      mcolon(pos(i), pos(i+1)-1-num+offset);
end
