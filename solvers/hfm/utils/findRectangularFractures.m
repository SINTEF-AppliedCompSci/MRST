function rectangular = findRectangularFractures(fractures)
% Indicates the rectangular fractures using a structure containing, for
% each fracture, the set of points defining the fracture polygon in 3D.

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

rectangular = [];
for i = 1:numel(fractures)
    if size(fractures(i).points,1) == 4
        p = [fractures(i).points;fractures(i).points(1,:)];
        d = sqrt(sum(diff(p,1).^2,2));
        if d(1) == d(3) && d(2) == d(4) && d(1) + d(2) == d(3) + d(4)
            rectangular = [rectangular;i]; %#ok
        end
    end
end
            