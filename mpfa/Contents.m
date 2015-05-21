% Routines supporting the MPFA-O method for the pressure equation.
%
% Files
%   computeMultiPointTrans.m   - Compute multi-point transmissibilities.
%   incompMPFA.m               - Solve incompressible flow problem (fluxes/pressures) using MPFA-O method.
%   matrixBlocksFromSparse.m   - Extract block-diagonal matrix elements from sparse matrix
%   examples/showGridEffects.m - Example demonstrating basic use of the MPFA-O pressure solver.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
