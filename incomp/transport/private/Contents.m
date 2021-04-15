% Files
%   assembleTransportSource    - Form source term vector from individual contributions
%   computeTransportSourceTerm - Compute source term contributions for transport
%   initFaceMob                - Initialize upwind face mobility saturation indices according to Darcy flux.
%   initTransport              - Compute input to transport solver from flow calculation results.
%   newtonRaphson2ph           - Solve non-linear equation F(s)=0 using Newton-Raphson method.
%   twophaseJacobian           - Residual and Jacobian of single point upwind solver for two-phase flow.
%   twophaseUpwBE              - Implicit single point upwind solver for two-phase flow, no gravity.
%   twophaseUpwBEGrav          - Implicit single point upwind solver for two-phase flow, including gravity.
%   twophaseUpwFE              - Explicit single point upwind solver for two-phase flow, no gravity.
%   twophaseUpwFEGrav          - Explicit single point upwind solver for two-phase flow, including gravity.

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
