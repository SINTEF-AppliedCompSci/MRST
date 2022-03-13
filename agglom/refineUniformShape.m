function p = refineUniformShape(p, G, indicator, NU, varargin)
%Refine blocks in a partition
%
% SYNOPSIS:
%   p = refineUniformShape(p, G, I, NU)
%   p = refineUniformShape(p, G, I, NU, 'cartDims', [n,m,..])
%   p = refineUniformShape(p, G, I, NU, 'fixPart', pf)
%
% DESCRIPTION:
%   Given a partition 'p' of a grid 'G', this function adds a fixed
%   subdivision of each coarse block in which the cummulative cell
%   indicators 'I' exceed an upper limit. Here, 'I' is interpreted as a
%   density and multiplied by the cell volume to derive total indicator
%   value for each block. No block should have a cummulative I-value larger
%   than NU/G.cells.num times the total cummulative I-value for the whole
%   grid.
%
% REQUIRED PARAMETERS:
%   p  - Partition vector
%   G  - Grid data structure discretising the reservoir model
%   I  - Cell-wise value of some measure/indicator function used for
%        deciding which blocks to refine.
%   NU - Upper bound on indicator values per block
%
% OPTIONAL PARAMETERS supplied in 'key'/value pairs ('pn1', pv1, ..):
%   cartDims - Rectangular subdivison based on the bounding box of each
%              coarse grid block
%   fixPart  - A predefined partition vector which is substituted into each
%              coarse block with a too large indicator. If the partition is
%              a matrix, the columns are applied in sequence until the
%              upper bound is fulfilled.
%   preserve - Boolean. If true, the original partition is preserved.
%
% RETURNS:
%   p - Partition vector after refining.

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

opt = struct('cartDims', repmat(2, [1, G.griddim]), ...
             'fixPart', [], ...
             'preserve', false);
opt = merge_options(opt, varargin{:});

assert(numel(opt.cartDims)==G.griddim, ...
      ['Number of Cartisan refinement factors must ', ...
       'match number of grid dimensions.\n', ...
       'Expected %d factors, but got %d.'], ...
       G.griddim, numel(opt.cartDims));

pCart = reshape(1:prod(opt.cartDims), opt.cartDims);

indicator      = indicator .* G.cells.volumes;
blockIndicator = accumarray(p, indicator);
upper_bound    = NU*sum(indicator)/G.cells.num;

% If a partition vector is provided, we first impose this partition vector
% in all blocks that exceed the upper limit. Likewise, if a partition array
% is provided, we apply the columns sequentially until the upper bound is
% fulfilled.
q = p;
if ~isempty(opt.fixPart)
   assert(any([numel(opt.fixPart),size(opt.fixPart,1)]==G.cells.num), ...
      'Wrong dimension of partition vector');
   if numel(opt.fixPart)==G.cells.num
      opt.fixPart = reshape(opt.fixPart, [], 1);
   end
   for i=1:size(opt.fixPart,2)
      bfix = false(1,max(p));
      bfix(blockIndicator>upper_bound) = true;
      p(bfix(p)) = opt.fixPart(bfix(p),i) + max(p);
      blockIndicator = accumarray(p, indicator);
   end
end
if opt.preserve,
    [~,~,p] = unique([p, q],'rows');
    p = compressPartition(p);
end

% Impose a rectangular partition using the bounding box for each
% block that exceeds the upper limit
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

   blks    = pick(p);
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
