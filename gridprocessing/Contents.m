% GRIDPROCESSING
%   Construct and manipulate the MRST grid datastructures.
%
% Files
%   buildCornerPtNodes      - Construct physical nodal coordinates for CP grid.
%   buildCornerPtPillars    - Construct physical nodal coordinates for CP grid.
%   cart2active             - Compute active cell numbers from linear Cartesian index.
%   cartGrid                - Construct 2d or 3d Cartesian grid in physical space.
%   cellNodes               - Extract local-to-global vertex numbering for grid cells.
%   checkAndRepairZCORN     - Detect and repair artifacts that may occur in corner-point specification.
%   computeGeometry         - Add geometry information (centroids, volumes, areas) to a grid
%   extended_grid_structure - Extended grid structure 
%   extractSubgrid          - Construct valid grid definition from subset of existing grid cells.
%   glue2DGrid              - Connect two 2D grids along common edges
%   grid_structure          - Grid structure used in the MATLAB Reservoir Simulation Toolbox.
%   grdeclXYZ               - Get corner-point pillars and coordinates in alternate format
%   hexahedralGrid          - Construct valid grid definition from points and list of hexahedra
%   makeInternalBoundary    - Make internal boundary in grid along specified faces.
%   makeLayeredGrid         - Extrude 2D Grid to layered 3D grid with specified layering structure
%   pebi                    - Compute dual grid of triangular grid G.
%   processFaults           - Construct fault structure from input specification (keyword `FAULTS`)
%   processGRDECL           - Compute grid topology and geometry from pillar grid description.
%   processNNC              - Establish explicit non-neighbouring connections from NNC keyword
%   processPINCH            - Establish vertical non-neighbouring across pinched-out layers
%   refineDeck              - Refine the grid resolution of a deck, and update other information
%   refineGrdecl            - Refine an Eclipse grid (`GRDECL file`) with a specified factor in each of
%   removeCells             - Remove cells from grid and renumber cells, faces and nodes.
%   removeFaultBdryFaces    - Remove fault faces on boundary
%   removeInternalBoundary  - Remove internal boundary in grid by merging faces in face list N
%   removeIntGrid           - Cast any grid fields that are presently int32 to double
%   removePinch             - Uniquify nodes, remove pinched faces and cells.
%   removeShortEdges        - Replace short edges in grid G by a single node.
%   splitDisconnectedGrid   - Split grid into disconnected components
%   tensorGrid              - Construct Cartesian grid with variable physical cell sizes.
%   tessellationGrid        - Construct valid grid definition from points and tessellation list
%   tetrahedralGrid         - Construct valid grid definition from points and tetrahedron list
%   triangleGrid            - Construct valid grid definition from points and triangle list
%   triangulateFaces        - Split face f in grid G into subfaces.

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
