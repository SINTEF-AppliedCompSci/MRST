% Support for wells and boundary conditions.
%
% Files
%   addBC               - Add boundary condition to (new or existing) BC object
%   addSource           - Add an explicit source to (new or existing) source object.
%   addWell             - Insert a well into the simulation model.
%   boundaryFaceIndices - Retrieve face indices belonging to subset of global outer faces.
%   controlWell         - Change well control target and/or value.
%   fluxside            - Impose flux boundary condition on global side.
%   pside               - Impose pressure boundary condition on global side.
%   psideh              - Impose hydrostatic Dirichlet boundary condition (pressure) on global side.
%   verticalWell        - Insert a vertical well into the simulation model.
%   processWells        - Construct MRST well structure from ECLIPSE input deck control.

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
