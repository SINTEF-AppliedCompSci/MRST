function [new_endp,ratio] = extend_frac2D(G, frac_endp,edges,out,cnum,cnodes,faces)
% extend_frac is used to extend a 2D fracture line when it only partially
% penetrates matrix cells at its end points. Following this, the ratio of
% true fracture length to its extended length inside the concerned matrix
% cell is returned. CI is computed for the extended fracture inside the
% corresponding matrix cell and scaled by this ratio. See Lee et al, Water
% Resources Research, 2001.

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


points = [G.cells.centroids(cnum,:);G.nodes.coords(cnodes,:);G.faces.centroids(faces,:)];
tri = delaunayTriangulation(points(:,1),points(:,2));
K = convexHull(tri);
edges2 = edges(~out.intAdjacencyMatrix,:); % Edges not intersected by the fracture
[xq,yq] = line_intersect(frac_endp,edges2);
[~,on] = inpolygon(xq,yq,tri.Points(K,1),tri.Points(K,2));
if all(on==0)
    cdist = sqrt((yq-G.cells.centroids(cnum,2)).^2 + (xq-G.cells.centroids(cnum,1)).^2);
    [~,on] = min(cdist);
end
new_endp = unique([xq(on),yq(on)],'rows');
frac_endp = [frac_endp(1:2:3).',frac_endp(2:2:4).'];
in = inpolygon(frac_endp(:,1),frac_endp(:,2),tri.Points(K,1),tri.Points(K,2));
true_endp = frac_endp(in,:);
startp = [uniquetol(out.intMatrixX(out.intAdjacencyMatrix),eps*100), ...
        uniquetol(out.intMatrixY(out.intAdjacencyMatrix),eps*100)];
if isempty(new_endp)
    error('End points not being recognized by inpolygon !');
end
new_endp = new_endp(~ismembertol(new_endp,startp,eps*100,'ByRows',true),:);
ratio = pdist_euclid([startp;true_endp])/pdist_euclid([startp;unique(new_endp,'rows')]);
%{
figure;triplot(tri); hold on; plotGrid(G,'FaceColor','none','EdgeAlpha',0.02); plot(frac_endp)
hold on; plot(frac_endp(:,1),frac_endp(:,2));
plot(xq,yq,'r*')
plot(tri.Points(K,1),tri.Points(K,2))
%}
return

function [xq,yq] = line_intersect(frac_endp,edges2)
frac_m = diff(frac_endp(2:2:4))/diff(frac_endp(1:2:3));
edge_m = zeros(size(edges2,1),1);
xq = zeros(size(edges2,1),1);
yq = zeros(size(edges2,1),1);
for i = 1:size(edges2,1)
    edge_m(i) = diff(edges2(i,2:2:4))/diff(edges2(i,1:2:3));
end
A = zeros(2,2); B = zeros(2,1);
if isinf(frac_m)
    A(1,1) = 1; A(1,2) = 0; B(1) = frac_endp(1);
    fractype = 1;
elseif frac_m == 0
    A(1,1) = 0; A(1,2) = 1; B(1) = frac_endp(2);
    fractype = 2;
else
    A(1,1) = -frac_m; A(1,2) = 1; B(1) = frac_endp(2)-frac_m*frac_endp(1); % y-m*x = y1-m*x1
    fractype = 3;
end
for i = 1:size(edges2,1)
    if isinf(edge_m(i))
        A(2,1) = 1; A(2,2) = 0; B(2) = edges2(i,1);
        edgetype = 1;
    elseif edge_m(i) == 0
        A(2,1) = 0; A(2,2) = 1; B(2) = edges2(i,2);
        edgetype = 2;
    else
        A(2,1) = -edge_m(i); A(2,2) = 1; B(2) = edges2(i,2)-edge_m(i)*edges2(i,1); % y-m*x = y1-m*x1
        edgetype = 3;
    end
    switch fractype
        case 1
            switch edgetype
                case 1
                    X = [NaN; NaN];
                case 2
                    X = [B(1); B(2)];
                otherwise
                    X = A\B;
            end
        case 2
            switch edgetype
                case 1
                    X = [B(2); B(1)];
                case 2
                    X = [NaN; NaN];
                otherwise
                    X = A\B;
            end
        otherwise
            X = A\B;
    end
    xq(i) = round(X(1),10);
    yq(i) = round(X(2),10);
end
