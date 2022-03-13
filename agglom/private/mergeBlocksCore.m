function mrg = mergeBlocksCore(bN, bIVol, bIFlw, bINum, lbnd, ubnd, NB, cfac)
%Core implementation of MERGE primitive
%
% SYNOPSIS:
%   mrg = mergeBlocksCore(bN, bIVol, bIFlw, bINum, lbnd, ubnd, NB, cfac)
%
% PARAMETERS:
%   bN    - Coarse-scale neighbourship definition, possibly derived using
%           function 'blockNeighbourship'.
%
%   bIVol - Accumulated volume indicator per block.  If a particular block
%           is empty, its 'bIVol' value should be NaN.
%
%   bIFlw - Relative flow indicator per block.  If a particular block is
%           empty, its 'bIFlw' value should be NaN.
%
%   bINum - Number of cells in each block
%
%   lbnd  - Lower block volume bound.  This algorithm will merge a block B
%           into a neighbouring block if bIVol(B)<lbnd.
%
%   ubnd  - Upper block flow bound.  The merging process will attempt to
%           honour the criterion bIFlw(B)*bIVol(B)<=ubnd.
%
%   NB    - Upper bound on number of cells within a single block.
%
%   cfac  - When violating the soft constraints on block flow (and number
%           of cells within a block) by a factor cfac, the constraints(s)
%           become hard constraint(s).
%
% RETURNS:
%   mrg   - A merging operator.  Specifically, a NUMEL(bIVol)-element
%           vector containing (a subset of) the numbers 1:NUMEL(bIVol) such
%           that
%
%               p2 = mrg(p)
%
%           in which 'p' is a partition vector will create a new, possibly
%           uncompressed, partition vector 'p2' defining the merged blocks.
%           Use function 'compressPartition' to remove empty blocks from
%           this block defintion.
%
% SEE ALSO:
%   `mergeBlocks2`, `blockNeighbourship`, `compressPartition`.

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


   assert(cfac>=1, 'mergeBlocksCore:AssertionFailed',...
      'Parameter <cfac> must be in the range [1,inf)');

   mrg = (1 : size(bIVol, 1)) .';

   conn = blockConnectivity(double(bN));

   bI = struct('Vol', bIVol, 'Flw', bIFlw, 'Num', bINum);

   for b = reshape(candidates(bI, lbnd), 1, []),
      if viable(b, conn, bI, lbnd),
         [mrg, conn, bI] = ...
            merge_into_neigh(b, mrg, conn, bI, ubnd, NB, cfac);
      end
   end
end

%--------------------------------------------------------------------------

function b = candidates(bI, lbnd)
   b      = find(bI.Vol < lbnd);
   [i, i] = sort(bI.Vol(b));              %#ok
   b      = b(i);
end

%--------------------------------------------------------------------------

function tf = viable(b, conn, bI, lbnd)
   i  = blockNeighbours(conn, b);
   tf = (bI.Vol(b) < lbnd) && any(isfinite(bI.Flw(i)));
end

%--------------------------------------------------------------------------

function [mrg, conn, bI] = ...
      merge_into_neigh(b, mrg, conn, bI, ubnd, NB, cfac)

   n   = blockNeighbours(conn, b);
   flw = bI.Vol(b)*bI.Flw(b) + bI.Vol(n).*bI.Flw(n);
   num = bI.Num(b) + bI.Num(n);

   feasible = ~( (flw > ubnd) | (num>NB) );

   if any(feasible),
      % Merge into feasible neighbour that most closely matches (block)
      % flow indicator of 'b'.

      n = n(feasible);
      [i, i] = min(abs(bI.Flw(b) - bI.Flw(n)));  %#ok
   else
      % Merge into neighbour that minimises violation of upper bounds.

      meas = flw/ubnd + num/NB;
      [i, i] = min(meas);                     %#ok
      if meas(i) > cfac, return, end;
   end

   into          = n(i);
   mrg(mrg == b) = into;

   % Update block indicator values for merging.
   % Recall:
   %   Vol and Num are additive while Flw is relative.
   bI.Flw(into) = bI.Vol(into)*bI.Flw(into) + ...
                  bI.Vol( b  )*bI.Flw( b  );

   bI.Vol(into) = bI.Vol(into) + bI.Vol(b)   ;   bI.Vol(b) = inf;
   bI.Flw(into) = bI.Flw(into) / bI.Vol(into);   bI.Flw(b) = inf;
   bI.Num(into) = bI.Num(into) + bI.Num(b)   ;   bI.Num(b) = inf;

   conn{into} = [conn{into}; conn{b}];
   conn{ b  } = [];
end
