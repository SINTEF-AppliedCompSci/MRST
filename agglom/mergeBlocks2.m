function p = mergeBlocks2(p, G, IVol, IFlw, NL, NU, varargin)
%Alternative implementation of Amalgamation 'MERGE' primitive
%
% SYNOPSIS:
%   p = mergeBlocks2(p, G, IVol, IFlw, NL, NU)
%   p = mergeBlocks2(p, G, IVol, IFlw, NL, NU, 'pn1', pv1, ...)
%
% PARAMETERS:
%   p      - Partition vector.  Possibly created by function
%            'segmentIndicator' or some other partitioning algorithm.  The
%            input partition vector should not contain any disconnected
%            blocks.  Function 'processPartition' will split such blocks.
%
%   G      - Grid structure.
%
%   IVol   - Block volume indicator.  One non-negative scalar, interpreted
%            as a density, for each cell in the grid 'G'.
%
%   IFlw   - Block flow indicator.  One non-negative scalar, interpreted as
%            a density, for each cell in the grid 'G'.  The algorithm will
%            merge a candidate block to its (feasible) neighbouring block
%            that most closely matches its own block flow.
%
%   NL, NU - Algorithm controlling parameters.  The algorithm will merge
%            blocks that violate the criterion
%
%                IVol(B) |B| >= (NL / n) IVol(G) |G|    (*)
%
%            while simultaneously attempting to uphold the upper bound flow
%            criterion
%
%                IFlw(B) |B| <= (NU / n) IFlw(G) |G|    (**)
%
%            Optionally, we may also try to uphold an upper bound on the
%            number of cells in each block
%
%                #cells(B) <= NB                        (***)
%
% OPTIONAL PARAMETERS:
%
%   nblock   - Upper bound on the number of cells in a single block.
%              Default: nblock = inf
%
%   cfac     - A relative factor at which the upper bound(s) are turned
%              into hard constraints. Default: cfac = inf
%
% RETURNS:
%   p      - Updated partition vector.  Typically contains fewer blocks
%            than the input partition vector.  None of the resulting blocks
%            should violate the criterion (*), but some of the blocks may
%            violate the criteria (**) and (***).
%
% SEE ALSO:
%   `mergeBlocks`, `refineBlocks`, `processPartition`.

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


   opt = struct('nblock', inf, 'cfac', inf);
   opt = merge_options(opt, varargin{:});

   bIVol = accumarray(p, IVol .* G.cells.volumes);
   bIFlw = accumarray(p, IFlw .* G.cells.volumes) ./ bIVol;
   bINum = accumarray(p, 1);

   lbnd = NL * sum(         bIVol) / G.cells.num;    % SUM(IVol .* volumes)
   ubnd = NU * sum(bIFlw .* bIVol) / G.cells.num;    % SUM(IFlw .* volumes)

   N    = getNeighbourship(G, 'Topological', false);
   bN   = blockNeighbourship(N, p);

   mrg  = mergeBlocksCore(bN, bIVol, bIFlw, bINum, ...
                          lbnd, ubnd, opt.nblock, opt.cfac);

   p = compressPartition(mrg(p));

   q = processPartition(G, p);
   if ~all(p==q),
      warning(msgid('Part:Disconnected'), ...
         'Resulting partition is disconnected. Repairing before continuing.');
      p = q;
   end

   if any(accumarray(p, IVol .* G.cells.volumes) < lbnd),
      warning(msgid('LBnd:Violated'), ...
         'Some blocks still violate lower (volume) bound after merging.');
   end
end
