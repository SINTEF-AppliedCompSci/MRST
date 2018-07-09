function d_avg = getAvgFracDist2D(G, fracp, cell_num,cnodes)
% getAvgFracDist computes the average normal distance of each point inside
% a matrix cell from every fracture line that intersects with this cell.
% This is the same as <d> used for computing conductivity index (CI, see
% CIcalculator2D) as described in SPE-103901-PA, Lyong Li and Seong H. Lee,
% 2008.

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


faces = G.cells.faces(G.cells.facePos(cell_num):G.cells.facePos(cell_num+1)-1,1);
t = 0:0.05:1;
xq = t*fracp(1,1)+(1-t)*fracp(2,1);
yq = t*fracp(1,2)+(1-t)*fracp(2,2);
xq = [xq';G.faces.centroids(faces,1);G.nodes.coords(cnodes,1)];
yq = [yq';G.faces.centroids(faces,2);G.nodes.coords(cnodes,2)];
pts = unique([xq,yq],'rows');
tri = delaunayTriangulation(pts(:,1),pts(:,2));
% triplot(tri,'k');

% Let fracture line be represented as Ax + By + C = 0
slope = diff(fracp(:,2))/diff(fracp(:,1));
if isinf(slope)
    A = 1; B = 0; C = -fracp(1,1);
elseif slope == 0
    A = 0; B = 1; C = -fracp(1,2);
else
    A = -slope;
    B = 1;
    C = -fracp(1,2)+fracp(1,1)*slope;
end
%
centroids = NaN(size(tri.ConnectivityList,1),2);
areas = NaN(size(tri.ConnectivityList,1),1);
d = NaN(size(tri.ConnectivityList,1),1);
dA = NaN(size(tri.ConnectivityList,1),1);
for i = 1:size(tri.ConnectivityList,1)
    pt_i = tri.Points(tri.ConnectivityList(i,:),:);
    centroids(i,:) = mean(pt_i);
    areas(i) = tri_area(pt_i(1,:),pt_i(2,:),pt_i(3,:));
    d(i) = abs(A*centroids(i,1)+B*centroids(i,2)+C)/(sqrt(A^2+B^2)); 
    dA(i) = d(i)*areas(i);
end
d_avg = sum(dA)/sum(areas);
return