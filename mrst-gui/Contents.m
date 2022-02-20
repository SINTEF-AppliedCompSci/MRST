% MRST-GUI
%
% Files
%   addFilters                 - Filter cells based on dataset.
%   colorizedHistogram         - Undocumented Utility Function
%   datasetSelector            - Create dataset selector either as standalone or integrated in gui.
%   editWells                  - Edit well setup
%   editStructGUI              - Undocumented Utility Function
%   extractSubcellsInteractive - Extract a subset of cells from a grid interactively.
%   fastRotateButton           - Hook in faster rotation for 3D plots.
%   getStructFields            - Dump possible plotting fields of a struct into human readable format.
%   interactiveSelection       - 
%   logColorbar                - Create colorbar for log dataset
%   mrstFigure                 - Create or select figure with convenient default options for 3d plotting
%   plotAdjustiblePlane        - {
%   plotCellVectorData         - Plot a vector field defined per cell
%   plotCurveTube              - Plot a curve as a tube with variable width
%   plotGridBarchart           - Plot barcharts on top of grid for cell data
%   plotToolbar                - Plot one or more datasets with interactive tools for visualization
%   plotWellData               - {
%   readStructField            - Read a field from a struct.
%   setRenderFixes             - Set fixes which may reduce crashes in MATLAB plots.
%   SummaryViewer              - Simple summary viewer. Usage (name without extension)
%   visualizeEclipseOutput     - Simple, experimental Eclipse output visualization function.

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
