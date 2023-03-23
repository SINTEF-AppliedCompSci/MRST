function area = tri_area(P1, P2, P3)
% tri_area(P1, P2, P3) calculates the triangle area given the coordinates
% of its vertices in P1, P2 and P3 using Heron's formula.
%
% Heron's Formula:
% s = semiperimeter
% A = sqrt(s * (s-a) * (s-b) * (s-c))
% Where a,b,c are lengths of the triangle edges

diff1 = P1 - P2;
diff2 = P1 - P3;
diff3 = P3 - P2;

a = norm(diff1);
b = norm(diff2);
c = norm(diff3);

% sort the elements
v = sort([a b c]);
a = v(3);
b = v(2);
c = v(1);

temp = b + c;
v1 = a + temp; % 2s
temp = a - b;
v2 = c - temp; % 2*(s-a)
v3 = c + temp; % 2*(s-b)
temp = b - c;
v4 = a + temp; % 2*(s-c)
area = 0.25 * sqrt(abs(v1*v2*v3*v4));

return

%{
Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.

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
