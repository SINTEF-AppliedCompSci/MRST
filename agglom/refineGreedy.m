function p = refineGreedy(p, G, IFlw, NU, varargin)
%Refine blocks in a partition using a greedy algorithm
%
% SYNOPSIS:
%   p = refineGreedy(p, G, IFlw, NU)
%   p = refineGreedy(p, G, IFlw, NU, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   The function performs a greedy refinement of all blocks in which the
%   (flow-based) indicator exceeds a prescribed upper bound. In each block,
%   the algorithm picks the cell that is furthest away from the block
%   center and starts growing a new block inward by adding all neighbouring
%   cells of the new block until the upper bound on the indicator is
%   exceeded. This process is repeated until the sum of the indicator
%   values of the remaining cells inside the original block is below the
%   threshold.
%
% PARAMETERS:
%   p      - Partition vector  Possibly created by function
%            'segmentIndicator' or some other partitioning algorithm.  The
%            input partition vector should not contain any disconnected
%            blocks.  Function 'processPartition' will split such blocks.
%
%   G      - Grid structure
%
%   IFlw   - Block flow indicator.  One non-negative scalar, interpreted as
%            a density, for each cell in the grid 'G'.  The algorithm will
%            merge a candidate block to its (feasible) neighbouring block
%            that most closely matches its own block flow.
%
%   NU     - Algorithm controlling parameter.  The algorithm will refine
%            blocks with too much flow--i.e., blocks for which
%
%                IFlw(B) |B| >= (NU / n) IFlw(G) |G|    (*)
%
% OPTIONAL PARAMETERS:
%
%   nlevel   - Specifies the definition of neighbourship used in the greedy
%              algorithm. Level 1 uses only the face neighbors, while level
%              2 also includes cells that have faces in common with at
%              least two of the level-1 neighbours.
%              Default: nlevel = 2
%
%   verbose  - Whether or not display number of blocks in the resulting
%              partition. Default value dependent upon global verbose
%              setting of function 'mrstVerbose'.
%
% RETURNS:
%   p      - Updated partition vector.  Typically contains more blocks
%            than the input partition vector.  Some of the resulting blocks
%            may violate the criterion (*) since the greedy algorithm adds
%            one ring of neighbors at the time when growing new blocks.
%
% SEE ALSO:
%   `refineGreedy2`, `refineGreedy3`, `refineBlocks`, `segmentIndicator`, `processPartition`

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

opt = struct('nlevel', 2, 'verbose', mrstVerbose);
opt = merge_options(opt, varargin{:});

% Construct neighbor matrix
NeighborMatrix = connectivity(G);

numBlocks = max(p); % counter for making new blocks

% Set indicator values, upper bound, etc
IFlw = IFlw.*G.cells.volumes;
ubnd = NU*sum(IFlw)/G.cells.num;
ui   = accumarray(p,IFlw);

%% Refining algorithm
%  Go through each coarse block in the partition and refine if necessary.
for i=reshape(find(ui > ubnd), 1, []),

   cellsInBlock = find(p==i);

   Ci = G.cells.centroids(cellsInBlock,:)';
   Ci = bsxfun(@minus, Ci, Ci(:,1));
   dc = sum(abs(Ci));  % distance vector of cell centers to block centers

   NM = NeighborMatrix(:,cellsInBlock);
   NM = NM(cellsInBlock,:);
   while any(dc)
      [maxVal, cellIdx] = max(dc);
      cellIdx           = min(cellIdx);

      [cIdxBlock, jj, cells] = find(cellsInBlock);

      e = zeros(length(cIdxBlock),1);
      e(cIdxBlock==cellIdx) = 1;

      I = NM(:,cIdxBlock);
      I = I(cIdxBlock,:);

      Ve = 0; Veo = -1;
      while (Ve > Veo) && (Ve < ubnd)
         e = I*e; % Find level-1 neighbors

         % Find level-2 neighbors
         if(opt.nlevel == 2)
            e        = e>0;
            e_neigh  = ((I*e)>0) - e;
            fe_neigh = find(e_neigh);

            no_neigh = I(fe_neigh,:)*e;
            e_neigh(fe_neigh(no_neigh==1))=0;

            e = e + e_neigh;  % Add the new neighbors to original "e"
         end

         Veo=Ve+1e-15;
         neighboringCells=find(e); ce=cells(neighboringCells);
         Ve=sum(IFlw(ce));
      end

      % Insert the new blocks at the end of the partition
      numBlocks=numBlocks+1;
      p(ce)= numBlocks;  % Add new block in the partition vector

      % Remove the processed cells from the block arrays
      cellsInNewBlock               = cIdxBlock(neighboringCells);
      dc(cellsInNewBlock)           = 0;
      cellsInBlock(cellsInNewBlock) = 0;
   end
end

 %% Update partition vector
 p = compressPartition(p);

end

%--------------------------------------------------------------------------

function A = connectivity(G)
   incBdry = false;  % Exclude boundary connections
   incDiag = true;   % Include diagonal entries in connectivity matrix
   N = getNeighbourship(G, 'Topological', incBdry);
   A = getConnectivityMatrix(N, incDiag, G.cells.num);
end
