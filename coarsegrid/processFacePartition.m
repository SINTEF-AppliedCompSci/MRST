function pf = processFacePartition(g, p, pf)
% Ensure that all coarse faces are connected collections of fine faces.
%
% SYNOPSIS:
%   pf = processFacepartition(g, p, pf)
%
% PARAMETERS:
%   g  - Grid structure as described by grid_structure.
%
%   p  - Partition vector of size [G.cells.num, 1] describing the coarse
%        grid.  We assume that all coarse blocks are connected.  The
%        partition vector is often created by function partitionCartGrid
%        or function partitionUI.
%
%   pf - Partition vector on faces of 'G'.
%
% RETURNS:
%   pf - Processed partition face vector.
%
% COMMENTS:
%   Currently only implemented for 2D grids.
%
% SEE ALSO:
%   `processPartition`, `cellPartitionToFacePartition`

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


   assert (isfield(g, 'griddim'), ...
          ['Function ''%s'' is only supported for grids which ', ...
           'declare their dimensions.'], mfilename);

   if g.griddim ~= 2,
      error(msgid('GridDim:Unsupported'), ...
           ['Function ''%s'' in only supported in ', ...
            'two space dimensions.'], mfilename);
   end

   % Identify faces on the boundary or between coarse blocks --------------
   P = [0; p];
   B = sort(P(g.faces.neighbors+1), 2);
   f = find(any(B == 0, 2) | (B(:,1) ~= B(:,2))) ;

   % face-nodes associated with f -----------------------------------------
   fnode_ix = mcolon(g.faces.nodePos(f), g.faces.nodePos(f + 1) - 1);

   % Create node -> face map for faces in list f --------------------------
   % For each node ni connected to ne of the faces in f, find all faces
   % in f also connected to ni.
   fno    = rldecode(f, g.faces.nodePos(f+1) - g.faces.nodePos(f));
   tmp    = sortrows([g.faces.nodes(fnode_ix), fno]);
   [N, n] = rlencode(tmp(:,1));                                        %#ok

   pos    = cumsum([1; n]);
   n2f    = tmp(:,2);
   clear tmp fno

   % Remove vertices with more than one partition tag ---------------------
   x   = [pf(n2f, :), B(n2f, :)];
   nno = rldecode(1:numel(n), n, 2)';
   i   = accumarray(nno, ...
            all(x(rldecode(pos(1:end-1), diff(pos)),:) == x, 2)) == n;
   n2f = n2f(rldecode(i, n));
   assert(all( n(i) == 2 ), ...
      'There should be no more than two faces for each connection...');

   % Assemble adjacency graph ---------------------------------------------
   m(f) = 1:numel(f);
   n2f  = double(reshape(n2f, 2, [])');
   a    = sparse(m(n2f), m(fliplr(n2f)), 1, numel(f), numel(f));
   clear m

   % Find connected coarse faces ------------------------------------------
   % Reassign face partition for faces in f, such that each coarse face is
   % a connected collection of fine grid faces.
   [p,q,r]  = dmperm(a+speye(size(a)));                                %#ok
   pf(f(p)) = rldecode(1:numel(diff(r)), diff(r), 2)';
end
