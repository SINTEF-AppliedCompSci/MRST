% INVENTORY
%
% Files
%   postprocessStates             - This function takes the results of a simulation ('states', a cell array of
%   postprocessStates3D           - This function is a wrapper that calls 'postprocessStates' on a result from a 3D
%   massTrappingDistributionVEADI - Compute the trapping status distribution of CO2 in each cell of a top-surface grid
%   plotTrappingDistribution      - Generate Trapping Inventory Plot From Simulation Result

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
