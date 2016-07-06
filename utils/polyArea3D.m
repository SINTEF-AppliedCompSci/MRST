function A = polyArea3D(p)
% Computes the surface area of a polygon defined by the set of points 'p'
% in 3D

assert(size(p,1)>=3,'3 or more points needed to define a plane!');
assert(all(iscoplanar(p)),'All Points passed must be coplanar');

pdash = [p;p(1,:)];
diffp = diff(pdash,1);
normal = cross(diffp(1,:), diffp(2,:));
normal = normal/norm(normal);
A = 0;
for i = 1:size(p,1);
    A = A + cross(pdash(i,:),pdash(i+1,:));
end
A = 0.5*dot(normal,A);
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