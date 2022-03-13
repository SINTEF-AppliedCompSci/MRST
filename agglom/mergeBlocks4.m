function p = mergeBlocks4(p, G, I, N, flist, varargin)
%Amalgamation 'MERGE' primitive adapted to fault information
%
% SYNOPSIS:
%   p = mergeBlocks4(p, G, I, N, flt)
%   p = mergeBlocks4(p, G, I, N, flt, 'pn1', pv1, ...)
%
% PARAMETERS:
%   p   - Partition vector.  Possibly created by function
%         'applySuccessivePart' or some other partitioning algorithm.  The
%         input partition vector should not contain any disconnected
%         blocks.  Function 'processPartition' will split such blocks.
%
%   G   - Grid structure.
%
%   I   - Indicator function structure. The following indicators
%         must be present as individual structure fields:
%
%            - volume  -- Block volume indicator, use unit value to get
%                         bulk volume and the porosity to get pore volume
%            - feature -- Block feature indicator (integer number)
%
%         The algorithm will generally merge a candidate block to the its
%         (feasible) neighbouring block that most closely matches its own
%         block flow.
%
%   N   - Algorithm controlling parameters.  Structure that features the
%         following fields
%
%            - volumeLow --
%                 Lower bound on total block volume.  The algorithm will
%                 merge blocks that violate the criterion
%
%                    I.volume(B) |B| >= (N.volumeLow / n) I.volume(G) |G|
%
%            - featHigh --
%                 Upper bound on total block flow.  While merging blocks
%                 that violate the lower block volume criterion, it will
%                 attempt to uphold the upper bound flow criterion
%
%                    I.flow(B) |B| <= (N.flowHigh / n) I.flow(G) |G|
%
%            - blocksHigh --
%                 Upper bound on total number of fine-scale cells in a
%                 block.  When merging blocks that violate the lower block
%                 volume criterion, the algorithm will attempt to uphold
%                 the criterion
%
%                    #cells(B) <= N.blocksHigh
%
%                 OPTIONAL.  Treated as INF (inactive) if not present.
%
%   flt - M-element structure array that defines faults.  Each structure
%         element must define the following fields:
%
%            - hard -- Numeric array of faces (face indices) that
%                 constitute the fault's connections.
%
%            - soft -- Numeric array of faces (face indices) that extend
%                 the fault towards the boundary.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   cfac - A relative factor at which the upper bound(s) are turned into
%          hard constraints.  Default: cfac = INF (do not turn criteria into
%          hard constraints).
%
%   merge_vert -
%          LOGICAL (Boolean) flag indicating whether or not to merge blocks
%          in the vertical direction.  Default value: merge_vert = TRUE (do
%          merge blocks in vertical direction).
%
% RETURNS:
%   p - Updated partition vector.  Typically contains fewer blocks than the
%       input partition vector.  None of the resulting blocks should
%       violate the criterion lower bound on block volumes, but some of the
%       blocks may violate the upper bounds on total block flow and/or
%       total number of cells in a block.
%
% SEE ALSO:
%   `mergeBlocks2`, `refineBlocks`, `processPartition`.

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

   opt = struct('cfac', inf, 'merge_vert', true);
   opt = merge_options(opt, varargin{:});

   bI    = block_indicators(G, p, I);
   bN    = define_neighbourship(G, p, flist, opt.merge_vert);
   bound = compute_bounds(G.cells.num, bI, N);

   %mrg = mergeBlocksCore(bN, bI, bound, opt.cfac);
   %mrg = (1 : size(bI.Vol, 1)) .';
   conn = blockConnectivity(double(bN));
   for b = reshape(candidates(bI, bound), 1, []),
      if bI.Vol(b) < bound.volumeLow,
         [p, conn, bI] = ...
            merge_into_neigh(p, G, I, b, conn, bI, bound, opt.cfac);
      end
   end

   p = compressPartition(p);
   q = processPartition(G, p);
   if ~all(p == q),
      warning(msgid('Part:Disconnected'), ...
             ['Resulting partition is disconnected. ', ...
              'Repairing before continuing.']);
      p = q;
   end

   if any(accumarray(p, I.volume .* G.cells.volumes) < bound.volumeLow),
      warning(msgid('LBnd:Violated'), ...
             ['Some blocks still violate lower (volume) ', ...
              'bound after merging.']);
   end
end

%--------------------------------------------------------------------------

function bI = block_indicators(G, p, I)
   feat  = accumarray(p, I.feature, [], @mode);
   accum = sparse(p, 1 : G.cells.num, 1) * ...
           [ (I.volume .* G.cells.volumes), ...
             double(I.feature ~= feat(p)), ...
             ones([G.cells.num, 1]) ];

   bI = struct('Vol',  accum(:,1), ...
               'Feat', feat, ...
               'Dev',  accum(:,2)./accum(:,3), ...
               'Num',  accum(:,3));
end

%--------------------------------------------------------------------------

function bN = define_neighbourship(G, p, flist, vertical)
   neigh_kind = 'Topological';  % Include explicit NNCs

   % Include boundary connections.  Needed for consistent indexing using
   % 'flist.hard'
   incBdry = true;

   N = getNeighbourship(G, neigh_kind, incBdry);

   pick = true([size(N, 1), 1]);

   if ~vertical,
      vert_conn_p      = false([max(G.cells.faces(:,2)), 1]);
      vert_conn_p(5:6) = true;

      is_vert = accumarray(G.cells.faces(:,1), ...
                           vert_conn_p(G.cells.faces(:,2)), ...
                           [size(N,1), 1], @any);

      pick = pick & ~is_vert;  clear is_vert vert_p
   end

   if ~isempty(flist),
      pick = pick & eliminate_fault_conn(N, p, flist);
   end

   bN = blockNeighbourship(N(pick, :), p);
end

%--------------------------------------------------------------------------

function bound = compute_bounds(nc, bI, N)
   avgvol = sum(bI.Vol) / nc;
   bound = struct( ...
      'volumeLow',   N.volumeLow * avgvol, ...
      'volumeAvg',   N.volumeAvg * avgvol, ...
      'devHigh',     N.featureDev, ...
      'cellHigh',    inf);

   if isfield(N, 'cellHigh'),
      bound.cellHigh = N.cellHigh;
   end
end

%--------------------------------------------------------------------------

function pick = eliminate_fault_conn(N, p, flist)
   % Eliminate hard fault connections as well as (soft) extensions of those
   % hard faults that happen to affect the blocks intersected by the hard
   % faults.  Extensions in blocks unaffected by hard faults are treated as
   % regular connections.

   p1 = [ 0 ; reshape(p, [], 1) ];
   nh = cellfun('prodofsize', { flist.hard });  % #hard fault conns
   ns = cellfun('prodofsize', { flist.soft });  % #soft fault conns

   % hard = [ block , (hard) fault ID ]
   %
   % Note "2 * nh" to account for block pairs.
   hard = [ p1(reshape(N(vertcat(flist.hard), :) .', [], 1) + 1), ...
            rldecode(1 : numel(flist), 2 * nh, 2) .' ];

   hard = unique(hard, 'rows'); % Expensive.  Hope 'hard' isn't too large.

   % soft = concatenation of soft extension faces, one copy for each row in
   % 'hard'.  Potentially large.
   soft = vertcat(flist(hard(:, 2)).soft);

   % i = 'soft' indices for which one of the connecting blocks is affected
   % by the corresponding hard fault.
   i = any(bsxfun(@eq, p1(N(soft, :) + 1), ...
                  rldecode(hard(:, 1), ns(hard(:, 2)))), 2);

   % Mark hard faults and block-restricted soft extensions for exclusion.
   pick = true([size(N, 1), 1]);
   pick([soft(i) ; vertcat(flist.hard)]) = false;
end

%--------------------------------------------------------------------------

function b = candidates(bI, bound)
   b      = find(bI.Vol < bound.volumeLow);
   [i, i] = sort(bI.Vol(b));              %#ok
   b      = b(i);
end

%--------------------------------------------------------------------------

function [p, conn, bI] = ...
   merge_into_neigh(p, G, I, b, conn, bI, bound, cfac)

   n   = blockNeighbours(conn, b);
   n = n(isfinite(bI.Feat(n)));
   if isempty(n),
      return
   end

   % Compute number of cells in block, block volume, block feature by
   % majority vote, and percentage of cells whose feature deviate from the
   % assigned block value. As a measure of convexity, we compute the ratio
   % of the change in area to the change in volume.
   num = bI.Num(b) + bI.Num(n);
   vol = bI.Vol(b) + bI.Vol(n);
   [feat, dev, dAdV] = deal(zeros(size(num)));
   cellsb = find(p==b);
   for i=1:numel(feat)
      cellsn  = find(p==n(i));    facesn = boundaryFaces(G, cellsn);
      cells   = [cellsb; cellsn]; faces  = boundaryFaces(G, cells);
      dAdV(i) = (sum(G.faces.areas(faces)) - sum(G.faces.areas(facesn))) ...
         / bI.Vol(n(i));
      feat(i) = mode(I.feature(cells));
      dev(i)  = sum( I.feature(cells)~= feat(i)) / num(i);
   end

   % Merge into neighbour that minimises violation of upper bounds.
   % However, do not merge if error measure exceeds hard constraint.
   meas   = max(vol - bound.volumeAvg, 0)/bound.volumeAvg + ...
            dev/bound.devHigh + ...
            max(num - bound.cellHigh,0)/bound.cellHigh;
   [j, j] = sort(dAdV);                                         %#ok<ASGLU>
   [i, i] = min(meas(j));                                       %#ok<ASGLU>
   if meas(j(i)) > cfac,
      return
   end
   into = n(j(i));
   p(p == b) = into;

   % Update block indicator values for merging.
   bI.Vol (into) = vol(i);   bI.Vol (b) = inf;
   bI.Feat(into) = feat(i);  bI.Feat(b) = inf;
   bI.Dev (into) = dev(i);   bI.Dev (b) = inf;
   bI.Num (into) = num(i);   bI.Num (b) = inf;

   conn{into} = [conn{into}; conn{b}];
   conn{ b  } = [];
end
