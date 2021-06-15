% UTILS
%
% Files
%   control2schedule         - Convert control vector u to schedule
%   evalObjective            - Objective (and gradient) evaluation function based on input control vector u
%   initSimpleScaledADIFluid - version of initSimpleScaledADIFluid with additional relperm scaling
%   scaleConstraints         - Linear constraint scaling
%   schedule2control         - Convert schedule to control vector
%   setupConstraints         - Setup linear constraints for scaled problem. Assumes linConst applies to

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
