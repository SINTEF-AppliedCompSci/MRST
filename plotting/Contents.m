% Routines for visual inspection of grid geometry and field properties.
%
% Files
%   boundaryFaces.m         - Extract boundary faces from set of grid cells.
%   outlineCoarseGrid.m     - Impose outline of coarse grid on existing grid plot.
%   plotBlockAndNeighbors.m - Plot a coarse block and its neighbors to current axes (reversed Z axis).
%   plotCellData.m          - Plot exterior grid faces, coloured by given data, to current axes.
%   plotContours.m          - Plot contours of cell data.
%   plotFaceData.m          - Plot face data on exterior grid faces to current axes (reversed Z axis).
%   plotFaces.m             - Plot selection of coloured grid faces to current axes (reversed Z axis).
%   plotFaults.m            - Plot faults in model
%   plotGrid.m              - Plot exterior grid faces to current axes (reversed Z axis).
%   plotGridVolumes.m       - Plot partially transparent isosurfaces for a set of values
%   plotSlice.m             - Plot Cartesian slices of cell data on faces
%   plotWell.m              - Plot well trajectories into current axes.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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
