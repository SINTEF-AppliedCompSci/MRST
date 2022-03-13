function p = refineRecursiveCart(p, G, indicator, NU, cartDims)
%Refine blocks by recursively applying Cartesian refinement pattern
%
% SYNOPSIS:
%   p = refineRecursiveCart(p, G, I, NU, cartDims)
%
% DESCRIPTION:
%   Given a partition 'p' of a grid 'G', this function recursively
%   subdivides each coarse block in which the accumulated block indicators
%   'I' exceed an upper limit.  The subdivision is done according to a
%   prescribed Cartesian refinement pattern.
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
%   NU - Upper bound on indicator values per block.
%
%   cartDims -
%        Rectangular subdivison based on the bounding box of each coarse
%        grid block.
%
% NOTE:
%   This function uses features of the 'gridtools' and 'coarsegrid' modules
%   so those modules must be active in order to use function
%   'refineRecursiveCart'.
%
% RETURNS:
%   p - Partition vector after refinement.  The partition vector contains
%       connected, non-empty blocks only because the final step of function
%       'refineRecursiveCart' is to pass the vector through the functions
%       'compressPartition' and 'processPartition' from the 'coarsegrid'
%       module.
%
% SEE ALSO:
%   `applySuccessivePart`, `sampleFromBox`, `compressPartition`, `processPartition`.

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

   assert (all(cartDims > 0), ...
           'All Cartesian refinement sizes must be strictly positive.');

   assert (numel(cartDims) == G.griddim, ...
          ['Number of Cartisan refinement factors must ', ...
           'match number of grid dimensions.\n', ...
           'Expected %d factors, but got %d.'], G.griddim, numel(cartDims));

   pCart = reshape(1 : prod(cartDims), cartDims);

   indicator   = indicator .* G.cells.volumes;
   upper_bound = NU * sum(indicator) / G.cells.num;

   % Impose a rectangular partition using the bounding box for each block
   % that exceeds the upper limit.
   mp      = max(p);
   violate = @(p) find(accumarray(p, indicator) > upper_bound);
   row     = @(v) reshape(v, 1, []);
   pick    = @(p) row(violate(p));
   blks    = pick(p);

   proceed = ~ isempty(blks);

   while proceed,
      split = false;

      for b = blks,
         cellsInBlock = find(p == b);

         part  = reshape(sampleFromBox(G, pCart, cellsInBlock), [], 1);
         n     = cumsum(accumarray(part, 1) > 0);
         split = split || (n(end) > 1);

         p(cellsInBlock) = part + mp;

         mp = mp + numel(n);
      end

      blks = pick(p);

      % Continue as long as there are any remaining blocks that violate the
      % upper bound *AND* we successfully split at least one block among
      % the previous block set.
      %
      proceed = split && ~ isempty(blks);
   end

   if ~ isempty(blks),
      warning(msgid('UBnd:Violated'), ...
             ['Final partition violates Upper Bound in %d blocks ', ...
              'after Cartesian Refinement.'], numel(blks));
   end

   % Compact the values and make sure that each block is singly connected
   p = processPartition(G, compressPartition(p));
end
