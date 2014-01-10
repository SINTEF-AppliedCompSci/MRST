function p = refineSuccessivePart(p, G, indicator, NU, pfixed)
%Refine blocks by successively applying static background partitions
%
% SYNOPSIS:
%   p = refineSuccessivePart(p, G, I, NU, pf)
%
% DESCRIPTION:
%   Given a partition 'p' of a grid 'G', this function applies a succession
%   of fixed subdivisions of each coarse block in which accumulated cell
%   indicators 'I' exceed an upper limit.
%
%   Here, 'I' is interpreted as a density and multiplied by the cell volume
%   to derive the total indicator value for each cell.  No block should
%   have a total accumulated indicator value exceeding "NU / G.cells.num"
%   times the total accumulated indicator value for the whole grid,
%
%       SUM(I .* G.cells.volumes)
%
% PARAMETERS:
%   p  - Partition vector.
%
%   G  - Grid data structure discretising the reservoir model.
%
%   I  - Cell-wise value of some measure/indicator function used for
%        deciding which blocks to refine.
%
%   NU - Upper bound on indicator values per block.  Positive scalar.
%
%   pf - A predefined partition vector which is substituted into each
%        coarse block with indicator values that are too large.  If the
%        fixed partition is a (G.cells.num)-by-m matrix, the columns are
%        each applied once, in sequence of increasing column numbers.
%
% NOTE:
%   This function uses features of the 'coarsegrid' module so that module
%   must be active in order to use function 'refineSuccessivePart'.
%
% RETURNS:
%   p - Partition vector after refinement.  The partition vector contains
%       connected, non-empty blocks only because the final step of function
%       'refineSuccessivePart' is to pass the vector through the functions
%       'compressPartition' and 'processPartition' from the 'coarsegrid'
%       module.
%
% SEE ALSO:
%   compressPartition, processPartition

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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

   if numel(pfixed) == G.cells.num,
      pfixed = reshape(pfixed, [], 1);
   end

   for i = 1 : size(pfixed, 2),
      bfix = false([1, max(p)]);
      bfix(blockIndicator > upper_bound) = true;

      p(bfix(p)) = pfixed(bfix(p), i) + max(p);

      blockIndicator = blkInd(p);
   end

   % Compact the values and make sure that each block is singly connected.
   p = processPartition(G, compressPartition(p));
end
