% UTILS
%
% Files
%   addLimiter            - Add limiter to ned or existing limiter object.
%   assignDofFromState    - Assign dofs from state (typically initial state). All dofs
%   dgBasis               - Undocumented Utility Function
%   getLimiter            - Get limiter for dG simulations. Currently only supports TVB and
%   getMinMax             - Undocumented Utility Function
%   plotDGBasis           - Undocumented Utility Function
%   polyDim               - Computes the dimension of the space of polynomials of degree k or less
%   setupOperatorsDG      - Undocumented Utility Function
%   velocityInterpolation - Construct velocity vector from face fluxes

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
