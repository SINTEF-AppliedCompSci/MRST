function p = refineGreedy3(p, G, IFlw, NU, varargin)
%Refine blocks in a partition using a greedy algorithm
%
% SYNOPSIS:
%   p = refineGreedy3(p, G, IFlw, NU)
%   p = refineGreedy3(p, G, IFlw, NU, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   The function performs a greedy refinement of all blocks in which the
%   (flow-based) indicator exceeds a prescribed upper bound. In each block,
%   the algorithm picks the cell that is furthest away from the block
%   center and starts growing a new block inward by adding all neighbours
%   (sorted by decreasing number of faces shared with cells in the block)
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
%            blocks with too much flow or too many cells--i.e., blocks
%            which meet either of the critera
%
%                IFlw(B) |B| >= (NU / n) IFlw(G) |G|    (*)
%                n_B         >= NU                      (**)
%
% OPTIONAL PARAMETERS:
%
%   nlevel   - Specifies the definition of neighbourship used in the greedy
%              algorithm. Level-2 neighbours are all cells that share at
%              least one node.
%              Default: nlevel = 2
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

opt = struct('nlevel', 2, 'verbose', mrstVerbose);
opt = merge_options(opt, varargin{:});

% Construct logical neighbor matrix
cellno = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
col    = 1 + (G.faces.neighbors(G.cells.faces(:,1), 1) == cellno);

cneigh = G.faces.neighbors(sub2ind([G.faces.num, 2], ...
                                   double(G.cells.faces(:,1)), col));

% Construct neighbor matrix
NeighborMatrix = connectivity(G, opt);

numBlocks = max(p); % counter for making new blocks

% Set indicator values, upper bound, etc
IFlw = IFlw .* G.cells.volumes;
ubnd = NU * sum(IFlw) / G.cells.num;
ui   = accumarray(p, IFlw);
un   = accumarray(p, 1);

% Refinement algorithm
%  Go through each coarse block in the partition and refine if necessary.
for i = reshape(find((ui > ubnd) | (un > NU)), 1, []),

   cellsInBlock = find(p == i);

   Ci = G.cells.centroids(cellsInBlock, :) .';
   Ci = bsxfun(@minus, Ci, Ci(:,1));
   dc = sum(abs(Ci));  % distance vector of cell centers to block centers

   ubound = ubnd;
   % ubound = ui(i) / ceil(ui(i)/ubnd); % try to equalize the split

   NM = NeighborMatrix(:, cellsInBlock);
   NM = NM(cellsInBlock,:);

   % Record local cell indices of new block.
   %
   % Note: This gets ovewritten/filled anew for each iteration/new coarse
   % block in the outer WHILE loop (while any(dc)).
   nCells = zeros(size(cellsInBlock));
   while any(dc)
      [cellIdx, cellIdx] = max(dc); %#ok
      cellIdx            = min(cellIdx);

      pick = false([G.cells.num + 1, 1]);

      % Local     ~   global
      [cIdxBlock, jj, cells] = find(cellsInBlock);  %#ok

      e = zeros([length(cIdxBlock), 1]);  % Discovered (BFS(1)) region

      nCells(1) = find(cIdxBlock == cellIdx);
      e(nCells(1)) = 1; % Start growing from 'cellIdx'.

      nnc = 1;  % Number of cells in new block.

      % Local BFS(1) discovery operator.
      I = NM(:, cIdxBlock);
      I = I(cIdxBlock, :);

      c = cells(nCells(1)); Ve = IFlw(c); pick(c + 1) = true;

      minIFlw = min(IFlw(cells));

      while ~(Ve + minIFlw > ubound) && (nnc < NU),
         % Include new cells in block according to the following selection:
         %   - Consider only the immediate neighbours of 'c' that have not
         %     already been included
         %   - For each of these neighbours, count the number of joint
         %     interfaces in common with the new coarse block that we're
         %     currently growing
         %   - Order the candidate cells, descendingly, according to the
         %     number of joint interfaces
         %   - Include new cells until neighbour list is exhausted or we
         %     exceed the 'ubound' on total block flow.

         % 1) Immediate neighbours of 'c'.
         eo    = e;
         e     = (I * e) > 0;
         index = find(e - eo);

         if isempty(index),
            % There are no neighbours that have not already been
            % discovered.  We have, essentially, exhausted the list of
            % candidate block cells and should just terminate the growing
            % process and proceed to the next coarse block that violates
            % the 'ubound'.
            break
         end

         % 2) Count number of joint interfaces.
         c  = cells(index);
         ix = mcolon(G.cells.facePos(  c  ), ...
                     G.cells.facePos(c + 1) - 1 );
         nf = accumarray(cellno(ix), double(pick(cneigh(ix) + 1)));

         if ~ any(nf > 0),
            % Candidate cells don't have *any* joint interfaces with cells
            % already included in new block.  Strange, but possible.
            % Terminate searching/growing, and go on to next refinement
            % step.
            break
         end

         % 3) Order candidate cells, descendingly, according to number of
         %    joint interfaces.
         [ii, ii]   = sort(nf(c), 'descend'); %#ok
         [c, index] = deal(c(ii), index(ii));

         % 4) Compute indicator for new cells and add from the start of the
         %    list until the cumulative indicator exceeds threshold or list
         %    of candidate cells is exhausted.
         cumInd = Ve + cumsum(IFlw(c));
         index  = index(~(cumInd > ubound));
         nAdd   = numel(index);

         if ~ nAdd,
            % We can include no other cells in this block because we've
            % already exceeded the ubound on total block flow.
            break
         end

         % Record the new cells (local cell numbers)
         nCells(nnc + 1 : nnc + nAdd) = index;

         % Update list of cells included in this new block to prepare for
         % defining new neighbours.
         pick(cells(index) + 1) = true;

         nnc = nnc + nAdd;  % Number of block cells.
         Ve = cumInd(nAdd); % Total block flow.
      end

      % Insert the new blocks at the end of the partition
      numBlocks = numBlocks + 1;

      p(cells(nCells(1 : nnc))) = numBlocks;

      % Remove the processed cells from the block arrays
      cellsInNewBlock               = cIdxBlock(nCells(1 : nnc));

      dc(cellsInNewBlock)           = 0;
      cellsInBlock(cellsInNewBlock) = 0;
   end
end

% Update partition vector
p = processPartition(G, compressPartition(p));

end

%--------------------------------------------------------------------------

function A = connectivity(G, opt)
   % Construct neighbor matrix
   if opt.nlevel == 2,
      N = neighboursByNodes(G);
   else
      incBdry = false;  % Exclude boundary connections
      N = getNeighbourship(G, 'Toplogical', incBdry);
   end

   incDiag = true; % Include diagonal entries of connectivity matrix.

   A = getConnectivityMatrix(N, incDiag, G.cells.num);
end
