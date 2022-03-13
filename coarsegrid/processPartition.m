function partition2 = processPartition(G, partition, facelist)
%Split disconnected coarse blocks into new blocks.
%
% SYNOPSIS:
%   q = processPartition(G, p)
%   q = processPartition(G, p, facelist)
%
% PARAMETERS:
%   G        - Grid structure as described by grid_structure.
%
%   p        - Vector of size G.cells.num-by-1 of initial cell-to-block
%              mappings.  Assumed to be a vector of positive integers.  The
%              coarse block numbers are preserved for all blocks that are
%              internally connected.  Consequently, if all coarse blocks
%              are connected then ALL(q == p).
%
%   facelist - List of faces across which the connections between cells
%              should be removed before the processing
%
% RETURNS:
%   q        - Updated partition with only connected coarse blocks.
%
% SEE ALSO:
%   `processFacePartition`, `partitionUI`, `partitionCartGrid`.

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


non_empty            = find(accumarray(partition, 1) > 0);
[b2cPos, b2c, locno] = invertPartition(partition);

maxBlk = non_empty(end);

% Find list of valid connections: first, disregard connections between
% cells that belong to different coarse blocks
p1 = [0; partition];
N  = p1(G.faces.neighbors + 1);
valid = N(:,1) == N(:,2);
clear N p1

% Then remove connections that have been specified explicitly by user
if nargin>=3
   if isnumeric(facelist)
      assert(min(facelist)>0 && max(facelist)<=G.faces.num, ...
         'Face number not within valid range');
      valid(facelist) = false;
   elseif islogical(facelist) && numel(facelist)==G.faces.num
      valid(facelist) = false;
   else
      error(msgid('FaceList:NonNumeric'), ...
         'Face list ''faces'' is not numeric or a boolean per face.')
   end
end

partition2 = repmat(-1, size(partition));
for k = 1 : numel(non_empty)
   b = non_empty(k);

   %  1) Identify cells and block-internal connections for block 'b'.
   c = b2c(b2cPos(b) : b2cPos(b+1) - 1);
   f = G.cells.faces(mcolon(G.cells.facePos( c ), ...
                            G.cells.facePos(c+1) - 1), 1);

   locN = G.faces.neighbors(f, :);
   locN = locno(locN(valid(f), :));

   %  2) Construct symmetric adjacency matrix for block 'b' (non-zero
   %     diagonal) connectivity.
   % Ref:
   %    <URL:http://blogs.mathworks.com/steve/2007/03/20/
   %         connected-component-labeling-part-3/#comments>,
   %    specifically comment No. 8.
   %
   i   = (1 : numel(c)) .';
   adj = sparse([locN(:,1); locN(:,2); i], ...
                [locN(:,2); locN(:,1); i],  1 );

   %  3) Compute the connected components in this coarse block by means
   %     of Dulmage-Mendelsohn permutation.
   %
   [p, r, r] = dmperm(adj);   nComp = numel(r) - 1;                    %#ok

   %  4) Assign new block numbers to the individual connected components
   %     found in step 3 while taking care to preserve the original block
   %     number of the first (and, possibly, only) block (i.e., component).
   %     This ensures that the output partition equals the input partition
   %     if no coarse blocks must be split.  As a special case, array
   %     'blkNo' must be empty when nbc==0 to avoid triggering an assertion
   %     in 'rldecode'.  This condition occurs if coarse block 'b' contains
   %     a single cell.
   %
   blkNo            = [b, maxBlk + (1 : nComp - 1)];
   partition2(c(p)) = rldecode(blkNo, reshape(diff(r), 1, []), 2);

   maxBlk = maxBlk + max(nComp - 1, 0);
end

assert (all(partition2 > 0));

end
