function p = applySuccessivePart(p, G, indicator, NU, pfixed)
%Refine blocks by successively applying static background partitions
%
% SYNOPSIS:
%   p = applySuccessivePart(p, G, I, NU, pf)
%
% DESCRIPTION:
%   Given a partition 'p' of a grid 'G', this function applies a set of
%   fixed subdivisions pf(:,1) to pf(:,end) successively to coarse blocks
%   in the grid. In each block, the next subdivision in the succession is
%   only applied if the accumulative indicator inside the block exceeds a
%   prescribed threshold. The indicator function 'I' is interpreted as a
%   density and multiplied by the cell volume to derive the total indicator
%   value for each cell. No block should have a total accumulated indicator
%   value exceeding "NU / G.cells.num" times the total accumulated
%   indicator value for the whole grid,
%
%       SUM(I .* G.cells.volumes)
%
% PARAMETERS:
%   p  - Partition vector.
%
%   G  - Grid data structure
%
%   I  - Cell-wise value of some indicator density function used to
%        decide whether the next partition should be applied inside
%
%   NU - Upper bound on indicator values per block.  Positive scalar.
%
%   pf - A predefined set of partition vectors which are substituted into
%        each coarse block with indicator values that are too large.  If
%        the fixed partition is a (G.cells.num)-by-m matrix, the columns
%        are each applied once, in sequence of increasing column numbers.
%
% NOTE:
%   This function uses features of the 'coarsegrid' module, which must be
%   loaded before using the function.
%
% RETURNS:
%   p -  Partition vector after refinement.  The partition vector contains
%        connected, non-empty blocks only because the final step of
%        function 'applySuccessivePart' is to pass the vector through the
%        functions 'compressPartition' and 'processPartition' from the
%        'coarsegrid' module.
%
% EXAMPLES:
%   G = computeGeometry(cartGrid([5 5]));
%   p1 = partitionCartGrid(G.cartDims, [2 1]);
%   p2 = partitionCartGrid(G.cartDims, [1 2]);
%   p = ones(G.cells.num,1);
%   p = applySuccessivePart(p, G, p, 12, [p1 p2]);
%   plotCellData(G,p,'EdgeColor','k');
%
% SEE ALSO:
%   `compressPartition`, `processPartition`

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


   require coarsegrid

   assert (isnumeric(NU) && (numel(NU) == 1) && (NU > 0), ...
           '%s: Upper bound ''NU'' must be a positive scalar', mfilename);

   indicator      = indicator .* G.cells.volumes;
   upper_bound    = NU * sum(indicator) / G.cells.num;
   blkInd         = @(p) accumarray(p, indicator);
   blockIndicator = blkInd(p);

   % If a partition vector is provided, we first impose this partition
   % vector in all blocks that exceed the upper limit.  Likewise, if a
   % partition array is provided, we apply the columns successively until
   % the upper bound is fulfilled.

   assert (any([numel(pfixed), size(pfixed, 1)] == G.cells.num), ...
          ['Fixed background partition array must contain one ', ...
           'value for each grid cell.']);

   if numel(pfixed) == G.cells.num
      pfixed = reshape(pfixed, [], 1);
   end

   for i = 1 : size(pfixed, 2)
      nblk = max(p);

      bfix = false([nblk, 1]);
      bfix(blockIndicator > upper_bound) = true;

      if ~any(bfix)
         % No blocks violate upper bound so there is no need to consider
         % further background partitions.  Terminate refinement to omit
         % moderatly costly and otherwise idempotent updates to partition
         % vector 'p' and derived 'blockIndicator'.
         break
      end

      pick    = bfix(p);
      p(pick) = p(pick) + nblk*pfixed(pick, i);

      p = compressPartition(p);
      blockIndicator = blkInd(p);
   end

   % Compact the values and make sure that each block is singly connected.
   p = processPartition(G, compressPartition(p));
end
