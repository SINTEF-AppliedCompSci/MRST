% PLOTTING
%   Routines for visual inspection of grid geometry and field properties.
%
% Files
%   boundaryFaces         - Extract boundary faces from set of grid cells.
%   colorbarHist          - Make colorbar with histogram on top
%   mrstColorbar          - Append a colorbar with an accompanying histogram to the current axis
%   outlineCoarseGrid     - Impose outline of coarse grid on existing grid plot.
%   plotBlockAndNeighbors - Plot a coarse block and its neighbours to current axes (reversed Z axis).
%   plotCellData          - Plot exterior grid faces, coloured by given data, to current axes.
%   plotContours          - Plot contours of cell data.
%   plotFaceData          - Plot face data on exterior grid faces to current axes (reversed Z axis).
%   plotFaces             - Plot selection of coloured grid faces to current axes (reversed Z axis).
%   plotFaults            - Plot faults in model
%   plotGrid              - Plot exterior faces of grid to current axes.
%   plotGridVolumes       - Plot partially transparent isosurfaces for a set of values
%   plotNodeData          - Plot data defined on nodes of grid
%   plotSlice             - Plot Cartesian slices of cell data on faces
%   plotWell              - Plot well trajectories into current axes.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
