% COARSEGRID
%
% Files
%   cellPartitionToFacePartition - Construct partition of all grid faces from cell partition.
%   compressPartition            - Renumber coarse block partitioning to remove any empty coarse blocks.
%   generateCoarseGrid           - Form coarse grid from partition of fine-scale grid.
%   partitionCartGrid            - Partition a Cartesian grid.
%   partitionLayers              - Partition grid uniformly in logical (I,J) direction, non-uniformly in K.
%   partitionMETIS               - Partition grid according to connection strengths
%   partitionUI                  - Partition grid uniformly in logical space.
%   partitionTensor              - Partition Logically Cartesian Grid Into Tensor Product Blocks
%   processFacePartition         - Ensure that all coarse faces are connected collections of fine faces.
%   processPartition             - Split disconnected coarse blocks into new blocks.
%   refineNearWell               - Partition a set of points based on proximity to a well point
%   sortSubEdges                 - Sorts subfaces in a 2D coarse grid.
%   subFaces                     - Extract fine-grid faces constituting individual coarse grid faces.

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
