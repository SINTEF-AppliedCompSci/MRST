% GRID_UTILS
%
% Files
%   createAugmentedGrid - Extend grid with mappings needed for the virtual element solver
%   createGridMappings  - Add preliminary mappings to be used in createAugmentedGrid
%   flipGrid            - Flip a grid (z->x, x->y, y->z)
%   full_grid_structure - Extended grid structure used in the vemmech module
%   padGrdecl           - Add padding to corner-point grid so it is embedded in a box
%   refineGrdeclLayers  - Refine a GRDECL structure in the vertical direction
%   sortEdges           - Sort edges in G.faces.edges counter-clockwise to face orientation
%   verticalGrdecl      - Transform GRDECL pillars into vertical pillars

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
