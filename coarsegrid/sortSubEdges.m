function cg = sortSubEdges(cg)
%Sorts subfaces in a 2D coarse grid.
%
% SYNOPSIS:
%   [cg, s] = sortSubEdges(cg)
%
% PARAMETERS:
%   cg  - Coarse grid.
%
% RETURNS:
%   cg  - Row permutation of e such that edge e(j(i),:) is followed by edge
%         e(j(i+1),:).
%
%   s   - Sign of each subface, i.e., whether or not the edges have been
%         flipped.
%
% EXAMPLE:
%   g = cartGrid([2,2]);
%   p = ones(4, 1);
%   cg = generateCoarseGrid(g, p);
%   cg.faces.fconn
%   cg = sortSubEdges(cg);
%   cg.faces.fconn

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

   assert( isfield(cg, 'parent'),        ...
      'SORTSUBFACES has meaning for coarse grids only!');
   assert(~isfield(cg.parent, 'parent'), ...
      'SORTSUBFACES does not support nested coarse grids!');
   assert(cg.griddim == 2,               ...
      'SORTSUBFACES has no meaning for 3D grids.');

   cg = sort_fconn(cg);
   cg = sort_coarse_block_faces(cg);
end

%--------------------------------------------------------------------------

function e = get_fconn_edges(cg)
% extract (g.faces.connPos(end)-1)-by-2 array of node numbers, each row are
% start and end nodes of a fine-grid face in the list cg.faces.fconn.

  ix = mcolon(cg.parent.faces.nodePos(cg.faces.fconn), ...
              cg.parent.faces.nodePos(cg.faces.fconn+1)-1);
  e  = reshape(cg.parent.faces.nodes(ix), 2, [])';
end

%--------------------------------------------------------------------------

function E = get_coarse_edges(cg)
% extract cg.faces.num-by-2 array of node numbers, each row are start and
% end nodes of a coarse face.

   e   = get_fconn_edges(cg);
   fno = rldecode(1:cg.faces.num, 2*diff(cg.faces.connPos), 2)';
   tmp = sortrows([fno, reshape(e', [], 1), (1:numel(fno))']);
   [nodes, num] = rlencode(tmp(:,1:2));
   E   = reshape(nodes(num==1, 2), 2, [])';
end

%--------------------------------------------------------------------------

function cg = sort_coarse_block_faces(cg)
   E  = get_coarse_edges(cg);
   [j, s] = sortedges(E(cg.cells.faces(:,1),:), cg.cells.facePos);
   assert(all(s));
   cg.cells.faces = cg.cells.faces(j,:);
end

%--------------------------------------------------------------------------

function cg = sort_fconn(cg)
% Sort fine faces in cg.faces.fconn

   % extract edge nodes from fine grid cg.parent
   e = get_fconn_edges(cg);

   [j, s, success] = sortedges(e, cg.faces.connPos);
   if any(~success),
      warning(['Some coarse faces (%d to be precise) seem to '...
               'be disconnected sets of fine faces.'], sum(~success));
   end

   % find coarse faces where cg.faces.neighbors is not equal
   % "cg.partition(g.faces.neighbors
   n   = cg.parent.faces.neighbors(cg.faces.fconn(j),:);
   n(s<0, :) = n(s<0, [2,1]);
   p   = [0; cg.partition];
   N   = p(n+1);
   fno = rldecode(1:numel(cg.faces.connPos)-1, diff(cg.faces.connPos), 2)';
   tmp = rlencode([N, fno]);
   S   = all(cg.faces.neighbors == tmp(:,1:2), 2);

   % flip sequence of fine faces
   j = flip (j, cg.faces.connPos, find(~S));

   % change sequence of fine-grid faces in definition of coarse face.
   cg.faces.fconn = cg.faces.fconn(j);
end

%--------------------------------------------------------------------------

function v = flip(v, pos, i)
% For each section i, reverse values in v(pos(i):pos(i+1)-1).

   if nargin == 1, i = 1:numel(pos)-1; end
   ix1 = mcolon(pos(i),     pos(i+1)-1,  1);
   ix2 = mcolon(pos(i+1)-1, pos(i),     -1);
   v(ix1) = v(ix2);
end
