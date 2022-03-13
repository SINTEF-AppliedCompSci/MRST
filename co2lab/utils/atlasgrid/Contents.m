% ATLASGRID
%
% Files
%   contourAtlas            - Plot contour lines in 3D for height data
%   convertAtlasTo3D        - Create GRDECL struct from CO2 storage atlas thickness/top data
%   convertAtlasToStruct    - Create GRDECL struct from thickness/top data from the CO2 Storage Atlas 
%   getAtlasGrid            - Get GRDECL grids and datasets for CO2 Atlas datasets
%   getBarentsSeaNames      - Returns the formation names present in the Barents Sea.
%   getNorthSeaNames        - Returns the formation names present in the North Sea.
%   getNorwegianSeaNames    - Returns the formation names present in the Norwegian Sea.
%   processAAIGrid          - Process aii grid meta data to a grid
%   readAAIGrid             - Read AIIGrid from file.
%   updateWithHeterogeneity - Update deck and petroinfo to include heterogeneous rock properties

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
