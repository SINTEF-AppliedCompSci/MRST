% Models for various simple, incompressible fluids.
%
% Files
%   initCoreyFluid.m       - Initialize incompressible two-phase fluid model (res. sat., an. rel-perm).
%   initSWOFFluid.m        - Initialize incompressible two-phase fluid model (res. sat., an. rel-perm).
%   initSatnumFluid.m      - Initialize incompressible two-phase fluid model with satnum and relperm
%   initSatnumRelPerm.m    - Construct two-phase relperm evaluation function.
%   initSimpleFluid.m      - Initialize incompressible two-phase fluid model (analytic rel-perm).
%   initSimpleFluidJfunc.m - Initialize incompressible two-phase fluid model (analytic rel-perm).
%   initSimpleFluidPc.m    - Initialize incompressible two-phase fluid model (analytic rel-perm).
%   initSingleFluid.m      - Initialize incompressible single phase fluid model.

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
