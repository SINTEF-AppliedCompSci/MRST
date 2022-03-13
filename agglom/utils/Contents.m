% Files
%   blockConnectivity            - Build block-to-neighbours map by transposing neighbourship definition
%   blockNeighbours              - Identify the neighbours of a particular coarse block.
%   blockNeighbourship           - Derive coarse-scale neighbourship from fine-scale information
%   coarse_sat                   - Converts a fine saturation field to a coarse saturation field, weighted
%   convertBC2Coarse             - Convert fine-scale boundary conditions to coarse scale.
%   convertRock2Coarse           - Create coarse-scale porosity field for solving transport equation.
%   convertSource2Coarse         - Accumulate fine-scale source terms to coarse scale
%   findConfinedBlocks           - Identify coarse blocks confined entirely within a single other block.
%   localTransmissibilityUpscale - Use local transmissibility upscaling to find non-negative transmissibilities
%   mergeBlocksByConnections     - Merge blocks based on connection strength
%   mergeSingleNeighbour         - Undocumented Internal Utility Function
%   processPartitionSoft         - Split disconnected coarse blocks into new blocks.
%   removeConfinedBlocks         - Remove singular confined blocks and expose groups of confined blocks
%   signOfFineFacesOnCoarseFaces - Identify fine-scale flux direction corresponding to coarse-scale outflow

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
