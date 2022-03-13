% UTILS
%
% Files
%   addFracFlowInvADIFluid       - Add a fluid function for computing the inverse of the water fractional
%   addPcOWInvADIFluid           - Add a fluid function for computing the inverse of the capillary pressure.
%   cartBlockMap                 - Find identical coarse blocks in a partition of a Cartesian grid. This
%   createBlockFluid             - Extracting the fluid for the current coarse cell only.
%   createFracFlowTablesFromDeck - Undocumented Utility Function
%   initADIFluidOW               - Make a structure representing an oil-water fluid. This is might be a
%   initADIFluidOWPolymer        - Make a structure representing a three-component fluid (water, oil,
%   makePeriodicCartesianGrid    - Undocumented Utility Function
%   simulateToSteadyStateADI     - Run simulation to a steady state solution is found using fully implicit
%   struct2args                  - Converts a structure to a cell array with both fieldnames and values,

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
