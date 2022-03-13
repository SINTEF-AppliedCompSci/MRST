function p = mergeBlocks(p, G, IVol, IFlw, NL, varargin)
% Merge blocks in a partitioning that are smaller than the given limit.
%
% SYNOPSIS:
%   p = mergeBlocks(p, G, IVol, IFlw, NL)
%   p = mergeBlocks(p, G, IVol, IFlw, NL, ...
%                           'static_partition', p2);
%   p = mergeBLocks(p, G, IVol, IFlw, NL, ...
%                           'verbose', true);
%
% DESCRIPTION:
%   This function merges too small blocks to a neigbhoring blocks. Which
%   blocks to merge are decided by "IVol" and INTO which blocks they
%   are merged, is decided by "IFlw".
%
% REQUIRED PARAMETERS:
%   p     - Partition vector
%
%   G     - Grid data structure discretising the reservoir model
%           (fine grid, geological model).
%
%   IVol  - Cell-wise value of some measure/indicator function used for
%           deciding which blocks to merge.
%
%   IFlw  - Cell-wise value of some measure/indicator function used for
%           deciding which neighboring block to merge into.
%
%   NL    - Algorithm controlling parameter. The algorithm will merge
%           blocks that violate the criterion
%
%                 IVol(B) |B| >= (NL / n) IVol(G) |G|    (*)
%
%           with the neighboring block that has the closest IFlw value.
%
% OPTIONAL PARAMETERS:
%
%   verbose - Whether or not display number of blocks in the resulting
%           partition. Default value dependent on the global verbose
%           settings of function 'mrstVerbose'.
%
%   static_partition - A partitioning that should be preserved throughout
%           the merging.
%
% RETURNS:
%   p     - Partition vector after merging.
%
% SEE ALSO:
%   `segmentPartition`, `refineBlocks`

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


opt = struct('verbose' , mrstVerbose, ...
             'static_partition', []);
opt = merge_options(opt, varargin{:});

% Just for safety, intersect p with the static partition
if(~isempty(opt.static_partition))
   [b,i,p]=unique([p, opt.static_partition],'rows');
   p = processPartition(G,p);
else
   assert (all(p == processPartition(G, p)), ...
      'Inconsistent input partition. Disconnected?');
end

% Compute indicators and bounds
IVol  = IVol.*G.cells.volumes;
IFlw  = IFlw.*G.cells.volumes;
lbnd = NL*sum(IVol)/G.cells.num;

% Block indicator for deciding which blocks to merge
blockIVol = accumarray(p, IVol);

% Block indicator for deciding which block to merge INTO
blockIFlw = accumarray(p, IFlw)./blockIVol;

cells_nothingToDo = zeros(length(p), 1);

%--------------------------------------------------------------------------
% Merging algorithm
%
% As long as we have blocks with too small indicator value, we attempt to
% find a neighbor with similar indicator value (may be different than the
% first indicator) and merge the block to the neighbor.

while( min(blockIVol) < lbnd)

   nb = max(p);
   p1 = [0; p];
   pN = p1(G.faces.neighbors+1);

   pN = pN(all(pN>0,2),:);         % Exclude the boundary face connections.
   pN = pN(pN(:,1) ~= pN(:,2), :); % Exclude block-internal connections.

   % Build (symmetric) connectivity matrix.
   EE = sparse([pN(:,1); pN(:,2); (1 : nb).'], ...
               [pN(:,2); pN(:,1); (1 : nb).'], 1);

   % Merge blocks that are too small
   for block = reshape(find(blockIVol<lbnd), 1, [])
      cellsInBlock = (p==block);

      ei        = EE(block,:); % find faces belonging to this block
      ei(block) = 0;           % remove itself
      neighb    = find(ei);    % find neighbors

      if(isempty(neighb))
         % If no neighbors to the block to be merge, nothing to do.
         cells_nothingToDo(cellsInBlock)=1;
         continue;
      end;

      % Find difference in blockIVol value for deciding which
      % neighbor to merge into. If more possible neighbors, take the one
      % with lowest cell number.
      di1 = abs(blockIFlw(neighb) - blockIFlw(block));
      mergeToBlock = min(neighb(abs(di1-min(di1)) < eps));

      cellsInNewBlock = (p == mergeToBlock);

      found=1;
      % Make sure we do not violate the static partition.
      if(~isempty(opt.static_partition))
         tag  = opt.static_partition(cellsInBlock);     tag  = tag(1);
         tag2 = opt.static_partition(cellsInNewBlock);  tag2 = tag2(1);

         if(tag ~= tag2)
            % If the two blocks do not have the same "tag", they are not
            % within the same block in the static partition. Loop over the
            % possible neighbors to find a block with the same "tag".
            found=0;
            neighb_2=neighb;
            for teller=1:nnz(neighb_2)-1
               neighb_2 = setxor(neighb_2, mergeToBlock);
               mergeToBlock = min(neighb_2);
               cellsInNewBlock = (p==mergeToBlock);
               tag2 = opt.static_partition(cellsInNewBlock); tag2 = tag2(1);

               if(tag2 == tag)
                  found=1;
                  break;
               end
            end
         end
      end
      if(~found)
         cells_nothingToDo(cellsInBlock) = 1;
         continue;
      end
      p(cellsInBlock) = mergeToBlock; % -> jump in p-vector

      EE(mergeToBlock, neighb) = EE(mergeToBlock, neighb) + EE(block, neighb);
      EE(neighb, mergeToBlock) = EE(neighb, mergeToBlock) + EE(neighb, block);
      EE(:, block) = 0;
   end

   p = compressPartition(p);

   IVol(cells_nothingToDo>0) = lbnd*1000; % "remove" cells from indicator array

   p = processPartition(G, p);
   blockIVol = accumarray(p, IVol);
   blockIFlw = accumarray(p, IFlw)./blockIVol;

end

dispif(opt.verbose, 'MergeBlocks: %d blocks\n', max(p));

end

