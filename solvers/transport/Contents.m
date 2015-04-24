% Routines for solving transport/saturation equation.
%
% Files
%   computeTransportSourceTerm - Compute source terms for transport
%   explicitTransport          - Explicit single point upwind transpot solver for two-phase flow.
%   explicitTransportNew       - Explicit single point upwind transport solver for two-phase flow.
%   implicitTransport          - Implicit single point upwind transport solver for two-phase flow.
%   implicitTransportCG        - Coarse grid implicit single point upwind transport for two-phase flow
%   implicitTransportNew       - Implicit single point upwind transport solver for two-phase flow.
%   implicitTransportReorder   -
%   initFaceMob                - Initialize upwind face mobility saturation indices according to Darcy flux.
%   initTransport              - Compute input to transport solver from flow calculation results.
%   newtonRaphson2ph           - Solve non-linear equation F(s)=0 using Newton-Raphson method.
%   newtonRaphsonJEA           - newtonRaphson -- Solves the water saturation equation with
%   thetaTransport             - Implicit single point upwind transport solver for two-phase flow.
%   twophaseJacobian           - Implicit single point upwind solver for two-phase flow, including gravity.
%   twophaseUpwBE              - Implicit single point upwind solver for two-phase flow, no gravity.
%   twophaseUpwBEGrav          - Implicit single point upwind solver for two-phase flow, including gravity.
%   twophaseUpwFE              - Explicit single point upwind solver for two-phase flow, no gravity.
%   twophaseUpwFEGrav          - Explicit single point upwind solver for two-phase flow, including gravity.
%   twophaseUpwReorder         - Single point upwind solver for Buckley-Leverett flow based on reordering.

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
