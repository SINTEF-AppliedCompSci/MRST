function [in,on] = inPolygon3D(q, p, varargin)
% Indicates if a set of query points 'q' in 3D lie inside or on a polygon
% defined by coplanar points 'p'.

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

% q = n x 3 array of query points
% p = coplanar points defining a polygon
opt = struct('tolerance', eps, 'normal', [], 'checkIfCoplanar', false);
opt = merge_options(opt, varargin{:});
tolerance = opt.tolerance;
normal = opt.normal;

in = zeros(size(q,1),1);
on = zeros(size(q,1),1);
Sign = zeros(size(q,1),size(p,1));
allp = [q;p];
if opt.checkIfCoplanar
    warning('off','MATLAB:triangulation:EmptyTri3DWarnId');
    assert(iscoplanar(allp(:,1),allp(:,2),allp(:,3)),...
        'All Points passed must be coplanar');
    warning('on','MATLAB:triangulation:EmptyTri3DWarnId');
end
points = [p;p(1,:)];
diffp = diff(points,1);
if isempty(normal)
    normal = cross(diffp(1,:), diffp(2,:));
    normal = normal/norm(normal); % Normal of plane represented by the polygon
end

for j = 1:size(p,1)
    boundNormal = cross(-diffp(j,:),normal);
    boundNormal = boundNormal/norm(boundNormal);
    
    A = boundNormal(1); B = boundNormal(2); C = boundNormal(3);
    D = -dot(boundNormal,mean(points(j:j+1,:),1));
    
    sign = A*q(:,1) + B*q(:,2) + C*q(:,3) + D;
    
    sign(abs(sign)>tolerance) = sign(abs(sign)>tolerance)./abs(sign(abs(sign)>tolerance));
    sign(isnan(sign)) = 0;
    sign(abs(sign)<eps*100) = 0;
    Sign(:,j) = sign;
end

for j = 1:size(q,1)
    sign = Sign(j,:);
    if sum(sign)==size(p,1)
        in(j) = 1;
    elseif any(sign==0)
        if all(sign(sign~=0)==1)
            in(j) = 1;
            on(j) = 1;
        end
    end
end
in = logical(in); on = logical(on);
return
