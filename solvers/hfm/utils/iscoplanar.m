function [isit, disit] = iscoplanar(pxyz,varargin)
% Checks if a set of points lie on the plane defined by points 'pxyz' and
% returns an indicator value for each point along with its normal distance
% from the plane.

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

opt = struct('normal', [], 'point', [], 'distTolerance', eps);
opt = merge_options(opt, varargin{:});
dtol = opt.distTolerance;
normal = opt.normal;
p = opt.point;
[x, y, z] = deal(pxyz(:,1),pxyz(:,2),pxyz(:,3));

% Test 1: Very Strict

warning('off','MATLAB:triangulation:EmptyTri3DWarnId');
tri = delaunayTriangulation(x,y,z);
if isempty(tri.ConnectivityList);
    isit = ones(size(x));
else
    pp = [x,y,z];
    isit = ~ismember(pp,tri.Points);
end
warning('on','MATLAB:triangulation:EmptyTri3DWarnId');

% Test 2: Distance Based
if isempty(p)
    p = [x(1),y(1),z(1)];
end
if isempty(normal)
    diffp = diff([x,y,z],1);
    normal = cross(diffp(1,:), diffp(2,:));
    normal = normal/norm(normal);
end

A = normal(1); B = normal(2); C = normal(3);
D = -dot(normal,p);

dist = abs((A*x + B*y + C*z + D)./norm(normal));
disit = dist<=dtol;
isit = logical(isit);
% disit(isit) = true;
return

