function checkIfCoplanar(fracplanes)
% Asserts that a given set of points (in 3D) are coplanar.

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

warning('off','MATLAB:triangulation:EmptyTri3DWarnId');
for i = 1:numel(fracplanes)
    x = [];
    y = [];
    z = [];
    for j = 1:size(fracplanes(i).points,1)
        x = [x;fracplanes(i).points(j,1)]; %#ok
        y = [y;fracplanes(i).points(j,2)]; %#ok
        z = [z;fracplanes(i).points(j,3)]; %#ok
    end
    assert(all(iscoplanar([x,y,z])),['Input Points must be coplanar. Point set ',...
        num2str(i),' is not!']);
end

warning('on','MATLAB:triangulation:EmptyTri3DWarnId');
return