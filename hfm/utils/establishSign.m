function Sign = establishSign(G, fracplanes, varargin)
% Assigns a value of 0 (if the node/edge lies on the plane), 1 (if the
% node/edge lies on the positive side of the plane) or -1 (if the node/edge
% lies  on the negative side of the plane) to every corner point and edge
% of the input polygon with respect to the normal vector of the polygon and
% the normal vector of its bounding planes.

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

opt = struct('tolerance', eps*100);
opt = merge_options(opt, varargin{:});
tol = opt.tolerance;

Sign = struct;
fcents = G.faces.centroids;
ncoords = G.nodes.coords;
for i = 1:numel(fracplanes)
    normal = fracplanes(i).normal;
    points = [fracplanes(i).points;fracplanes(i).points(1,:)];
    p = mean(points,1);
    A = normal(1); B = normal(2); C = normal(3);
    D = -dot(normal,p);
    
    sign = A*fcents(:,1) + B*fcents(:,2) + C*fcents(:,3) + D;
    sign(abs(sign)>tol) = sign(abs(sign)>tol)./abs(sign(abs(sign)>tol));
    sign(isnan(sign)) = 0;
    sign(abs(sign)<tol) = 0;
    Sign(i).FaceSign(:,1) = sign;
    
    sign = A*ncoords(:,1) + B*ncoords(:,2) + C*ncoords(:,3) + D;
    sign(abs(sign)>tol) = sign(abs(sign)>tol)./abs(sign(abs(sign)>tol));
    sign(isnan(sign)) = 0;
    sign(abs(sign)<tol) = 0;
    Sign(i).NodeSign(:,1) = sign;
    
    for j = 1:size(fracplanes(i).boundNormals,1)
        normal = fracplanes(i).boundNormals(j,:);
        A = normal(1); B = normal(2); C = normal(3);
        D = -dot(normal,mean(points(j:j+1,:),1));
        
        sign = A*fcents(:,1) + B*fcents(:,2) + C*fcents(:,3) + D;
        sign(abs(sign)>tol) = sign(abs(sign)>tol)./abs(sign(abs(sign)>tol));
        sign(isnan(sign)) = 0;
        sign(abs(sign)<tol) = 0;
        Sign(i).FaceSign(:,j+1) = sign;
        
        sign = A*ncoords(:,1) + B*ncoords(:,2) + C*ncoords(:,3) + D;
        sign(abs(sign)>tol) = sign(abs(sign)>tol)./abs(sign(abs(sign)>tol));
        sign(isnan(sign)) = 0;
        sign(abs(sign)<tol) = 0;
        Sign(i).NodeSign(:,j+1) = sign;
    end
end

return