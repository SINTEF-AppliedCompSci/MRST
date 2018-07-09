function fracplanes = getPlaneNormals(fracplanes,varargin)
% Computes the normal vector for a set of coplanar points in 3D.
% Additionally, it returns the normal vector of the bounding planes
% perpendicular to the edges of the polygon defined by the set of input
% points.

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

opt = struct('normalize',false);
opt = merge_options(opt, varargin{:});

for i = 1:numel(fracplanes)
    points = [fracplanes(i).points;fracplanes(i).points(1,:)];
    
    diffp = diff(points,1);
    normal = cross(diffp(1,:), diffp(2,:));
    if opt.normalize
        normal = normal/norm(normal);
    end
    
    nBoundPlanes = size(points,1)-1;
    boundNormals = zeros(nBoundPlanes,3);
    for j = 1:nBoundPlanes
        boundNormals(j,:) = cross(-diffp(j,:),normal);
        if opt.normalize
            boundNormals(j,:) = boundNormals(j,:)/norm(boundNormals(j,:));
        end
    end
    fracplanes(i).normal = normal;
    fracplanes(i).boundNormals = boundNormals;
end