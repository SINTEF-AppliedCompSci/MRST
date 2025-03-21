% HELPERS
%
% Files
%   compareWellrates         - Undocumented Utility Function
%   getVEColors              - Return colors used for plotting VE models
%   make_testgrid            - Make a simple test grid, sloping along the x-direction and with a slight y
%   poroFromPerm             - Compute porosity from permeability using the inverse Cozeny-Karman
%   selectedResultsMultiplot - background variables can be: h|h_max|pressure|rs|sGmax|s|totalCO2
%   sloping_aquifer          - Undocumented Utility Function

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
