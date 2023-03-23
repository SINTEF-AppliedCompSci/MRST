% TESTGRIDS
%   Simple to moderately advanced synthetic test models
%
% Files
%   cartesianGrdecl      - Construct Cartesian grid with variable physical cell sizes.
%   createPinchedColumn  - Create a single column containing a single pinched layer of thickess 0.01
%   extrudedTriangleGrid - Build a synthetic grid with a curved fault in the lateral direction
%   makeModel3           - Build a synthetic geometry with two faults.
%   oneSlopingFault      - Make a GRDECL structure for a box grid with a single sloping fault.
%   pinchedLayersGrdecl  - Make a GRDECL structure for simple corner-point grid, possibly faulted.
%   pinchedNode          - Define two-cell corner-point specification with single, pinched vertex
%   pinchMiddleCell      - Create corner-point descriptions with variable number of pinched nodes
%   raisedColumn         - Create corner-point description of 2-by-1-by-2 grid with one fault
%   simpleGrdecl         - Make a GRDECL structure for simple corner-point grid, possibly faulted.
%   threeLayers          - Construct a corner point discretization of a three-layered structure.
%   twister              - Permutes x- and y-coordinates of nodes in a grid.

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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
