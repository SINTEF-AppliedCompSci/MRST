function [blks, pn] = findConfinedBlocks(G, p)
%Identify coarse blocks confined entirely within a single other block.
%
% SYNOPSIS:
%   blks       = findConfinedBlocks(G, p)
%   [blks, pn] = findConfinedBlocks(G, p)
%
% PARAMETERS:
%   G  - Grid structure.
%
%   p  - Partition vector defining coarse blocks.
%
% RETURNS:
%   blks - List of blocks that are confined inside another block, i.e.,
%          have only one neigbour in the interior of the grid.
%          NB! In the case of recursively confined blocks, only the
%          outermost and the innermost blocks will be listed in blks.
%
%   pn   - New partition in which each confined block is merged with the
%          block that surrounds it. NB! Will not work for nested
%          confinement.
%
% EXAMPLE:
%   % single confined block
%   G = cartGrid([4 4]);
%   p = ones(4,4); p(2,2) = 2; p(3:4,4)=3; p(4,2:3)=4;
%   blks = findConfinedBlocks(G,p(:));
%
%   % recursively confined blocks
%   G = cartGrid([7 4]);
%   p = ones(7,4); p(2:6,1:3)=2; p(3:5,1:2) = 3; p(4,1)=4;
%   [b, pn] = findConfinedBlocks(G,p(:));
%
% SEE ALSO:
%   `processPartition`, `blockNeighbourship`

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


   pn = p;
   bN = blockNeighbourship(G.faces.neighbors, p);
   b  = bN(:);
   a  = accumarray(b,1);
   bno = unique(b);
   blks = bno(a==1);

   if isempty(blks) || nargout==1
      return,
   elseif nargout<1
      warning(msgid('Part:Confined'), ...
         ['Detected blocks that are confined within other blocks; ', ...
         'please proceed with caution.']);
      return,
   end

   % Repair
   while ~isempty(blks)
      ix = false(max(b),1);
      ix(blks) = true;
      ix = reshape(ix(bN), [], 2);
      j = any(ix,2);
      bnew = sum(bN(j,:) .* ~ix(j,:), 2);
      [buf, i] = sort(sum(bN(j,:) .* ix(j,:), 2));
      bno(a==1) = bnew(i);
      pn = compressPartition(bno(pn)+1);       % +1 to avoid zeros in bno

      bN = blockNeighbourship(G.faces.neighbors, pn);
      b  = bN(:);
      a  = accumarray(b,ones(numel(b),1));
      bno = unique(b);
      blks = bno(a==1);
   end
end
