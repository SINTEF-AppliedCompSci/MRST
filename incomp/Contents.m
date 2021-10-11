% Support for incompressible flow/transport problems
%
% Files
%   capPressureRHS                - Compute capillary pressure contribution to system RHS
%   checkDrivingForcesIncomp      - Check if a solution is necessary for incompressible solvers. 
%   computeFacePressure           - Compute face pressure using two-point flux approximation.
%   computeIncompWellPressureDrop - Compute incompressible connection pressure drop for a single well
%   computePressureRHS            - Compute right-hand side contributions to pressure linear system.
%   incompTPFA                    - Solve incompressible flow problem (fluxes/pressures) using TPFA method.
%   treatLegacyForceOptions       - Internal function for ensuring that both W and Wells are supported as

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
