function p = psia()
%Compute numerical value, in units of Pascal, of one Psi.
%
% SYNOPSIS:
%   p = psia()
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   p - Numerical value, in units of Pascal (Pa), of a pressure of one Psi.
%
% NOTE:
%   The primary purpose of this utility function is to make examples
%   easier to read.

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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


   p = lbf() / inch()^2;
end
