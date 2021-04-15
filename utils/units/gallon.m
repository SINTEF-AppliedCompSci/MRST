function ga = gallon()
%Compute numerical value, in units of m^3, of one U.S. liquid gallon.
%
% SYNOPSIS:
%   ga = gallon()
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   ga - Numerical value, in units of m^3, of a volume of one U.S. liquid
%        gallon (== 231 cubic inches).
%
% NOTE:
%   The U.S. liquid gallon was historically defined as the volume of a
%   straight, vertical cylinder with an inner diameter of 7 inches and a
%   height of 6 inches.  This is 231 inch^3 if we assume that pi = 22/7.

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

   ga = 231 * inch^3;
end
