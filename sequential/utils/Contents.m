% UTILS
%
% Files
%   computeSequentialFluxes                 - Undocumented Utility Function
%   convertPressureBoundaryConditionsToFlux - Undocumented Utility Function
%   getPressureTransportIterations          - Get pressure, transport and outer loop iterations for sequential or
%   getSaturationUpwind                     - Get upwind flags for viscous and gravity parts of flux function.
%   getSequentialModelFromFI                - For a given fully implicit model, output the corresponding pressure/transport model
%   hybridUpwind                            - Hybrid upwinding - seperate gravity/capillary and viscous part
%   multiphaseUpwindIndices                 - Implementation of algorithm from "UPSTREAM DIFFERENCING FOR

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
