function p = refineGreedy4(p, G, IFlw, NU, varargin)
%Refine blocks in a partition using a greedy algorithm
%
% SYNOPSIS:
%   p = refineGreedy4(p, G, IFlw, NU)
%   p = refineGreedy4(p, G, IFlw, NU, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   The function performs a greedy refinement of all blocks in which the
%   (flow-based) indicator exceeds a prescribed upper bound. In each block,
%   the algorithm picks the cell that is furthest away from the block
%   center and starts growing a new block inward by adding all neighbours
%   (sorted by ascending difference of IFlw from the block average)
%   until the upper bound on the indicator is exceeded. This process is
%   repeated until the sum of the indicator values of the remaining
%   cells inside the block is below the threshold.
%
%   NB! The method uses 'neighboursByNodes' and may therefore be *slow*
%   e.g., compared to 'refineGreedy' and 'refineGreedy2'.
%
% PARAMETERS:
%   p      - Partition vector.  Possibly created by function
%            'segmentIndicator' or some other partitioning algorithm.  The
%            input partition vector should not contain any disconnected
%            blocks.  Function 'processPartition' will split such blocks.
%
%   G      - Grid structure.
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
%              algorithm. Level-2 neighbours are all cells that share at
%              least one node.
%              Default: nlevel = 2
%
%   nadd     - Upper bound on the number of cells to add during each
%              growing cycle. Default: nadd = inf
%
%   verbose  - Whether or not display number of blocks in the resulting
%              partition. Default value dependent upon global verbose
%              setting of function 'mrstVerbose'.
% RETURNS:
%   p      - Updated partition vector.  Typically contains more blocks
%            than the input partition vector.  None of the resulting blocks
%            should violate the criterion (*).
%
% SEE ALSO:
%   `refineGreedy`, `refineGreedy2`, `mergeBlocks`, `refineBlocks`, `processPartition`.

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

opt = struct('nlevel', 2, 'verbose', mrstVerbose, 'nadd', inf);
opt = merge_options(opt, varargin{:});

% Construct neighbor matrix
NeighborMatrix = connectivity(G, opt);

numBlocks = max(p); % counter for making new blocks

% Set indicator values, upper bound, etc
IFlw = IFlw.*G.cells.volumes;
ubnd = NU*sum(IFlw)/G.cells.num;
ui   = accumarray(p,IFlw);

% Refining algorithm
%  Go through each coarse block in the partition and refine if necessary.
for i=reshape(find(ui > ubnd), 1, []),

   cellsInBlock = find(p==i);

   Ci = G.cells.centroids(cellsInBlock,:)';
   Ci = bsxfun(@minus, Ci, Ci(:,1));
   dc = sum(abs(Ci));  % distance vector of cell centers to block centers

   ubound = ubnd;
   % ubound = ui(i) / ceil(ui(i)/ubnd); % try to equalize the split

   NM = NeighborMatrix(:,cellsInBlock);
   NM = NM(cellsInBlock,:);

   nCells = zeros(size(cellsInBlock));
   while any(dc)
      [maxVal, cellIdx] = max(dc);
      cellIdx           = min(cellIdx);
      pick = false([G.cells.num+1,1]);

      [cIdxBlock, jj, cells] = find(cellsInBlock);
      e = zeros(length(cIdxBlock),1);
      nCells(1) = find(cIdxBlock == cellIdx);
      e(nCells(1)) = 1;
      nnc = 1;

      I = NM(:,cIdxBlock);
      I = I(cIdxBlock,:);

      c = cells(nCells(1)); Ve = IFlw(c); pick(c+1) = true;
      minIFlw = min(IFlw(cells));
      eo = e;
      while ~(Ve + minIFlw > ubound)
         % Find neighbours, N(e)
         e  = (I*eo)>0; en = e - eo;
         index = find(en);
         if isempty(index), break, end

         % Compute indicator for new cells and add from the start of the
         % list until the cummulative indicator exceeds threshold. Sorting
         % candidate cells according to decreasing number of joint faces.
         c = cells(index);
         myInd = abs(Ve/nnc - IFlw(c));
         [nn,ii]=sort(myInd,'ascend');
         %if ~nn(1), disp('tilslag'), break, end
         index = index(ii); index = index(1:min(end,opt.nadd));
         c = c(ii);         c = c(1:min(end,opt.nadd));

         cumInd = Ve + cumsum(IFlw(c));
         index = index(~(cumInd > ubound));
         nAdd = numel(index);
         if ~nAdd, break, end
         nCells(nnc+1:nnc+nAdd) = index;
         pick(cells(index)+1) = true;
         nnc = nnc+nAdd;
         Ve = cumInd(nAdd);
         eo(nCells(1:nnc))=1;
      end

      % Insert the new blocks at the end of the partition
      numBlocks=numBlocks+1;
      p(cells(nCells(1:nnc))) = numBlocks;

      % Remove the processed cells from the block arrays
      cellsInNewBlock               = cIdxBlock(nCells(1:nnc));
      dc(cellsInNewBlock)           = 0;
      cellsInBlock(cellsInNewBlock) = 0;
   end
end

% Update partition vector
p = processPartition(G, compressPartition(p));

end

%--------------------------------------------------------------------------

function A = connectivity(G, opt)
   if opt.nlevel == 2,
      N = neighboursByNodes(G);
   else
      incBdry = false;  % Exclude boundary connections
      N = getNeighbourship(G, 'Toplogical', incBdry);
   end

   incDiag = true; % Include diagonal entries of connectivity matrix.

   A = getConnectivityMatrix(N, incDiag, G.cells.num);
end
