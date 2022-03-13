% Files
%   callMetisMatrix  - Partition connectivity graph whilst accounting for connection strengths
%   coarseDataToFine - Convert coarse grid dataset into fine grid representation
%   coarsenBC        - Construct coarse-grid boundary conditions from fine-grid boundary cond.
%   coarsenFlux      - Compute net flux on coarse faces
%   coarsenGeometry  - Add geometry (centroids, face normals, areas, ...) to a coarse grid
%   fineToCoarseSign - Compute sign change between fine faces and coarse faces.
%   invertPartition  - Invert partition (cell->block mapping) to create block->cell mapping.
%   isSamePartition  - Check if two partition vectors represent the same grid partition

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
