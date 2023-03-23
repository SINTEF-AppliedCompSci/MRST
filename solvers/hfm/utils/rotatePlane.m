function p = rotatePlane(points, normal, varargin)
% Rotates a set of points (assumed to be coplanar) in order to align with a
% plane defined by the input normal vector. Rotation is along the axis
% defined by the intersection of the two planes. If the plane defined by
% the input set of points is parallel to the input normal vector, then this
% function returns the points without making a transformation.

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

opt = struct('useNormal',[]);
opt = merge_options(opt,varargin{:});
assert(size(points,1)>=3,'3 or more points needed to define a plane!');
% assert(all(iscoplanar(points)),'All Points passed must be coplanar');
if isempty(opt.useNormal)
    diffp = diff(points,1);
    normal2 = cross(diffp(1,:), diffp(2,:));
    normal2 = normal2/norm(normal2); 
else
    normal2 = opt.useNormal/norm(opt.useNormal);
end
normal = normal/norm(normal);
if isequal(abs(normal),abs(normal2))
    p = points;
    return
end
rotAngle = -acos(dot(normal2,normal)); %rotation angle
rotAxis = -cross(normal,normal2); %rotation axis

M = makehgtform('axisrotate', rotAxis, rotAngle);
p = points*M(1:3,1:3);

return