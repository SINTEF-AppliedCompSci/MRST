function p = refineUniform(p, G, indicator, NU, varargin)
%Refine blocks in a partition by uniform partitioning
%
% SYNOPSIS:
%   p = refineUniform(p, G, indicator, NU)
%   p = refineUniform(p, G, indicator, NU, ...
%                            'cartDims',[nx_c ny_c nz_c])
%
% DESCRIPTION:
%   This function refines too large blocks by subdividing the blocks
%   uniformly according to the dimensions given in 'cartDims' (default
%   value is [2,2,2]). Calls partitionUI for the subdivision.
%
% REQUIRED PARAMETERS:
%   p         - Partition vector
%
%   G         - Grid data structure discretising the reservoir model
%               (fine grid, geological model).
%
%   indicator - Cell-wise value of some measure/indicator function used for
%               deciding which blocks to refine.
%
%   NU       - Upper bound
%
% OPTIONAL PARAMETERS:
%
%   cartDims - Dimensions of subdivison of the blocks.
%
%   verbose  - Whether or not display number of blocks in the resulting
%              partition. Default value dependent upon global verbose
%              settings of function 'mrstVerbose'.
%
% RETURNS:
%   p - Partition vector after refining.

% SEE ALSO:

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

if G.griddim==2,
   opt = struct('cartDims', [2, 2], 'verbose',  mrstVerbose);
else
   opt = struct('cartDims', [2, 2, 2], 'verbose',  mrstVerbose);
end
opt = merge_options(opt, varargin{:});

assert(isfield(G, 'cartDims'));

if isfield(G.cells, 'ijkMap'),
   ijk = G.cells.ijkMap;
else
   [ijk{1:G.griddim}] = ind2sub(G.cartDims, G.cells.indexMap(:));
   ijk = [ijk{:}];
end

indicator      = indicator .* G.cells.volumes;
blockIndicator = accumarray(p, indicator);
upper_bound    = NU*sum(indicator)/G.cells.num;

% While there are too large blocks, then refine
while max(blockIndicator) > upper_bound,
   [v, block]   = max(blockIndicator);
   cellsInBlock = find(p==block);

   % Build local ijk-map
   ijk_local = ijk(cellsInBlock,:);
   ijk_local = bsxfun(@minus, ijk_local, min(ijk_local)) + 1;

   % Generate a local grid for the current block, such that we can apply
   % partitionUI on it.
   G_loc = extractSubgrid(G, cellsInBlock);
   G_loc.cells.ijkMap = ijk_local;
   G_loc.cartDims     = max(ijk_local);

   if G_loc.griddim == 2,
      indexMap=sub2ind(G_loc.cartDims, ijk_local(:,1), ijk_local(:,2));
   else
      indexMap=sub2ind(G_loc.cartDims, ...
         ijk_local(:,1), ijk_local(:,2), ijk_local(:,3));
   end

   G_loc.cells.indexMap = indexMap;
   % G_loc = computeGeometry(G_loc);

   part = partitionUI(G_loc, min([opt.cartDims; G_loc.cartDims]));
   part = part + max(p);
   p(cellsInBlock) = part;

   % Updating indicator value for the blocks
   blockIndicator = accumarray(p, indicator);
end

% Update the partition vector
p = processPartition(G, compressPartition(p));

end
