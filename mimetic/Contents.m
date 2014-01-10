% Routines supporting the mimetic method for the pressure equation.
%
% Files
%   assembleWellSystem         - Generate pressure linear system components for wells.
%   calculateFaultTrans        - Compute transmissibilities of faults and add them to grid and optionaly
%   computeMimeticIP           - Compute mimetic inner product matrices.
%   packageWellSol             - Convert well fluxes and pressures to structure form.
%   solveIncompFlow            - Solve incompressible flow problem (fluxes/pressures).
%   solveIncompFlowFault       - Solve incompressible flow problem (fluxes/pressures) with fault and
%   unpackWellSystemComponents - Extract hybrid linear system components from wells.

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
