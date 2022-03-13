function pn = removeConfinedBlocks(G, p)
%Remove singular confined blocks and expose groups of confined blocks
%
% SYNOPSIS:
%   pn = removeConfinedBlocks(G, p)
%
% PARAMETERS:
%   G  - Grid structure.
%
%   p  - Partition vector defining coarse blocks.
%
% RETURNS:
%   pn - New partition in which each single confined block is merged with
%        the block that surrounds it and each block that completely
%        surrounds a set of other blocks is split in two (along a plane
%        orthogonal to the x-axis through the block center).  The resulting
%        partition is not necessarily singly connected and may have to be
%        processed by a call to function processPartition.
%
% EXAMPLE:
%   % Compare partitions before/after confined block removal.
%   % Intended for visual inspection only.
%   comparePartitions = @(p0,p1) ...
%      [rot90(p0) , zeros([size(p0,2),1]) , rot90(reshape(p1,size(p0)))]
%
%   % 1) Merge a singular confined block into the surrounding block.
%   G = cartGrid([4 4]);
%   p = ones(4,4); p(2,2) = 2; p(3:4,4)=3; p(4,2:3)=4;
%   pn = removeConfinedBlocks(G,p(:));
%   comparePartitions(p,pn)
%
%   % 2) Merge a set of individual blocks that are recursively confined
%   %    into the outermost surrounding block.
%   G = computeGeometry(cartGrid([7 7]));
%   p = ones(7,7); p(2:6,2:6)=2; p(3:5,3:5)=3;p(4,4)=4;
%   pn = removeConfinedBlocks(G,p(:));
%   comparePartitions(p,pn)
%
%   % 3) Expose a group of confined blocks.
%   G = computeGeometry(cartGrid([7 7]));
%   p = ones(7,7);p(:,4:7)=2;p(2:6,2:6)=3;p(3:5,3)=4;p(3:5,4)=5;p(3:5,5)=6;
%   pn = removeConfinedBlocks(G,p(:));
%   comparePartitions(p,pn)
%
% NOTE:
%   Function 'removeConfinedBlocks' uses the 'biconnected_components'
%   function of the MATLAB BGL library.  Consequently, if the MATLAB BGL is
%   not available then function 'removeConfinedBlocks' will fail.
%
% SEE ALSO:
%   `biconnected_components`, `processPartition`.

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


   assert (exist('biconnected_components', 'file') == 2, ...
           'The MATLAB BGL must be available to use function ''%s''.', ...
           mfilename);

   pn = single_confined(G, p);
   a  = group_confined(G, pn);

   if ~isempty(a)
      % Repair - in lack of a general solution: just split each block once
      % using a plane orthogonal to the x-axis through the centroid of the
      % coarse block
      mp = max(pn);
      for i=1:numel(a)
         cellsInBlock = find(pn==a(i)-1);
         x = G.cells.centroids(cellsInBlock,1);
         ce = cellsInBlock(x>mean(x));

         pn(ce)= mp+1;
         mp = mp+1;
      end
   end
end

%--------------------------------------------------------------------------

% Detect sets of isloated, single confined blocks
function p = single_confined(G, p)
   [bN, blks, bno, a] = single_neighbour_blocks(G, p);

   % Repair
   while ~isempty(blks)
      ix = false(max(max(bN)),1);
      ix(blks) = true;
      ix = reshape(ix(bN), [], 2);
      j = any(ix,2);
      bnew = sum(bN(j,:) .* ~ix(j,:), 2);
      [buf, i] = sort(sum(bN(j,:) .* ix(j,:), 2)); %#ok<ASGLU>
      bno(a==1) = bnew(i);
      p = compressPartition(bno(p)+1);       % +1 to avoid zeros in bno

      [bN, blks, bno, a] = single_neighbour_blocks(G, p);
   end
end

%--------------------------------------------------------------------------

% Detect groups of confined blocks
function a = group_confined(G, pn)
   bN = blockNeighbourship(G.faces.neighbors, pn, true)+1;
   C = sparse(bN, fliplr(bN), 1);
   a = biconnected_components(C);
end

%--------------------------------------------------------------------------

function [bN, blks, bno, a] = single_neighbour_blocks(G, p)
   bN = blockNeighbourship(G.faces.neighbors, p);
   b  = bN(:);
   a  = accumarray(b,1);
   bno = unique(b);
   blks = bno(a==1);
end
