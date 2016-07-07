% GRIDTOOLS
%   Functions that collect common operations on the MRST grid structure
%
% Files
%   checkGrid             - UNTITLED Summary of this function goes here
%   compareGrids          - Determine if two grid structures are the same.
%   connectedCells        - Compute connected components of grid cell subsets.
%   findEnclosingCell     - Find cells with closest centroid (in Euclidian norm) in a 2D grid.
%   getCellNoFaces        - Get a list over all half faces, accounting for possible NNC
%   getConnectivityMatrix - Derive global, undirected connectivity matrix from neighbourship relation.
%   getNeighbourship      - Retrieve neighbourship relation ("graph") from grid
%   gridAddHelpers        - Add helpers to existing grid structure for cleaner code structure.
%   gridCellFaces         - Find faces corresponding to a set of cells
%   gridCellNo            - Construct map from half-faces to cells or cell subset
%   gridCellNodes         - Extract nodes per cell in a particular set of cells
%   gridFaceNodes         - Find nodes corresponding to a set of faces
%   gridLogicalIndices    - Given grid G and optional subset of cells, find logical indices.
%   indirectionSub        - Index map of the type G.cells.facePos, G.faces.nodePos...
%   makePlanarGrid        - Construct 2D surface grid from faces of 3D grid.
%   neighboursByNodes     - Derive neighbourship from common node (vertex) relationship
%   removeNodes           - {
%   sampleFromBox         - Sample from data on a uniform Cartesian grid that covers the bounding box
%   sortGrid              - Permute nodes, faces and cells to sorted form
%   transform3Dto2Dgrid   - Transforms a 3D grid into a 2D grid.
%   translateGrid         - Move all grid coordinates according to particular translation
%   volumeByGaussGreens   - Compute cell volume by means of Gauss-Greens' formula

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
