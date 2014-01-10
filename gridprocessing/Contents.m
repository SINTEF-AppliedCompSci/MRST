% Examples demonstrating the construction and manipulation of grid datastructures.
%
% Files
%   buildCornerPtNodes.m     - Construct physical nodal coordinates for CP grid.
%   buildCornerPtPillars.m   - Construct physical nodal coordinates for CP grid.
%   cart2active.m            - Compute active cell numbers from linear Cartesian index.
%   cartGrid.m               - Construct 2d or 3d Cartesian grid in physical space.
%   cellNodes.m              - Extract local-to-global vertex numbering for grid cells.
%   computeGeometry.m        - Compute geometry of grid.
%   extractSubgrid.m         - Construct valid grid definition from subset of existing grid cells.
%   grid_structure.m         - Grid structure used in MATLAB Reservoir Simulation Toolbox.
%   hexahedralGrid.m         - Construct valid grid definition from points and list of hexahedra
%   makeInternalBoundary.m   - Make internal boundary in grid along FACES
%   makeLayeredGrid.m        - Extrude 2D grid to layered 3D grid with n layers.
%   pebi.m                   - Compute dual grid of triangular grid G.
%   processFaults.m          - Construct fault structure from input specification (keyword 'FAULTS')
%   processGRDECL.m          - Compute grid topology and geometry from pillar grid description.
%   removeCells.m            - Remove cells from grid and renumber cells, faces and nodes.
%   removeFaultBdryFaces.m   - Remove fault faces on boundary
%   removeInternalBoundary.m - Remove internal boundary in grid by merging faces in face list N
%   removePinch.m            - Uniquify nodes, remove pinched faces and cells.
%   tensorGrid.m             - Construct Cartesian grid with variable physical cell sizes.
%   tetrahedralGrid.m        - Construct valid grid definition from points and tetrahedron list
%   triangleGrid.m           - Construct valid grid definition from points and triangle list

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
