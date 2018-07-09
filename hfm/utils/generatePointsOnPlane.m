function [p,in,on] = generatePointsOnPlane(points,varargin)
% Generates points on a 3D convex polygon defined by a coplanar set of
% points

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

opt = struct('tolerance', 0, 'normal', [], 'ratio', 1e-2);
opt = merge_options(opt, varargin{:});
tolerance = opt.tolerance;
normal = opt.normal;

assert(size(points,1)>=3,'3 or more points needed to define a plane!');
assert(all(iscoplanar(points)),'All Points passed must be coplanar');

if isempty(normal)
    diffp = diff([points;points(1,:)],1);
    normal = cross(diffp(1,:), diffp(2,:)); % Not normalized
end

A = normal(1); B = normal(2); C = normal(3);
D = -dot(normal,points(1,:));
flag = 1*(A==0) + 2*(B==0) + 4*(C==0);

range = max(points)-min(points);
[~,sortInd] = sort(range);
if flag == 0
    flag = sortInd(1);
elseif flag == 1
    sortInd = sortInd(sortInd~=1); % Don't divide by A
    flag = sortInd(1);
elseif flag == 2
    sortInd = sortInd(sortInd~=2); % Don't divide by B
    flag = sortInd(1);
elseif flag == 4
    sortInd = sortInd(sortInd~=3); % Don't divide by C
    flag = sortInd(1);
elseif flag == 3 % // to x-y plane
elseif flag == 6 % // to y-z plane
elseif flag == 5 % // to x-z plane
end
t = transpose(0:opt.ratio:1);
switch flag
    case {1,6}
        yy = t*min(points(:,2)) + (1-t)*max(points(:,2));
        zz = t*min(points(:,3)) + (1-t)*max(points(:,3));
        [Y,Z] = meshgrid(yy,zz);
        X = (B*Y + C*Z + D)/(-A);
    case {2,5}
        xx = t*min(points(:,1)) + (1-t)*max(points(:,1));
        zz = t*min(points(:,3)) + (1-t)*max(points(:,3));
        [X,Z] = meshgrid(xx,zz);
        Y = (A*X + C*Z + D)/(-B);
    case {3,4}
        xx = t*min(points(:,1)) + (1-t)*max(points(:,1));
        yy = t*min(points(:,2)) + (1-t)*max(points(:,2));
        [X,Y] = meshgrid(xx,yy);
        Z = (A*X + B*Y + D)/(-C);
end

% Points on boundary
pdash = [points;points(1,:)];
addp = [];
for i = 1:size(points,1)
    addp = [addp;(1-t).*pdash(i,1) + t.*pdash(i+1,1),...
                 (1-t).*pdash(i,2) + t.*pdash(i+1,2),...
                 (1-t).*pdash(i,3) + t.*pdash(i+1,3)]; %#ok
end
query = unique([X(:),Y(:),Z(:);addp],'rows');

% Projecting on plane, not necessary
basisCoefficients= bsxfun(@minus,query,pdash(1,:))*null(normal);
q = bsxfun(@plus,basisCoefficients*null(normal).', pdash(1,:));

[init,on] = inPolygon3D(q, points, 'tolerance', tolerance, 'normal', normal);
p = query(init,:);
in = init(init);
on = on(init);

return
