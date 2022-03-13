% NUC Support for non-uniform coarsening/agglomeration method
%
% Files
%   applySuccessivePart - Refine blocks by successively applying static background partitions
%   mergeBlocks         - Merge blocks in a partitioning that are smaller than the given limit.
%   mergeBlocks2        - Alternative implementation of Amalgamation 'MERGE' primitive
%   mergeBlocks3        - Amalgamation 'MERGE' primitive adapted to fault information
%   mergeBlocks4        - Amalgamation 'MERGE' primitive adapted to fault information
%   refineBlocks        - Refine blocks for which indicator value exceeds given limit
%   refineGreedy        - Refine blocks in a partition using a greedy algorithm
%   refineGreedy2       - Refine blocks in a partition using a greedy algorithm
%   refineGreedy3       - Refine blocks in a partition using a greedy algorithm
%   refineGreedy4       - Refine blocks in a partition using a greedy algorithm
%   refineRecursiveCart - Refine blocks by recursively applying Cartesian refinement pattern
%   refineUniform       - Refine blocks in a partition by uniform partitioning
%   refineUniformShape  - Refine blocks in a partition
%   segmentIndicator    - Segments a fine grid into blocks according to indicator.

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
