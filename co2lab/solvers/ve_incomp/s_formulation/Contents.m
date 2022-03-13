% S_FORMULATION
%
% Files
%   gravPressureVE_s         - Computes innerproduct of (face_centroid - cell_centroid) * g for each face
%   initSimpleVEFluid_s      - Initialize incompressible two-phase fluid model for vertical average
%   primitivesMimeticVE_s    - Internal helper for topSurfaceGrid. Used to override mimetic primitives
%   twophaseJacobianWithVE_s - Residual and Jacobian of single point upwind solver for two-phase flow.

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
