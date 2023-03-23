%% Graphical View of the Grid Structure
% In this example, we will show details of the grid structures for three
% different grids: a Cartesian grid with one cell removed, a triangular
% grid, and a Voronoi grid. The examples will show cell numbers, node
% numbers, and face numbers.

%% REGULAR GRID
G = removeCells( cartGrid([3,2]), 2);
G = computeGeometry(G);

%% Plot cell, face, and node numbers
newplot;
plotGrid(G,'FaceColor',[0.95 0.95 0.95]); axis off;
hold on;
text(G.cells.centroids(:,1)-0.04, ...
   G.cells.centroids(:,2), num2str((1:G.cells.num)'),'FontSize',20);
plot(G.cells.centroids(:,1),...
   G.cells.centroids(:,2),'ok','MarkerSize',24);
text(G.faces.centroids(:,1)-0.045, G.faces.centroids(:,2), ...
   num2str((1:G.faces.num)'),'FontSize',16);
text(G.nodes.coords(:,1)-0.075, ...
   G.nodes.coords(:,2), num2str((1:G.nodes.num)'),'FontSize',18);
plot(G.nodes.coords(:,1),...
   G.nodes.coords(:,2),'sk','MarkerSize',24);
hold off;

%% Output content of cells.faces, faces.nodes, and faces.neighbors to file
faces =[ rldecode(1 : G.cells.num,diff(G.cells.facePos), 2).' G.cells.faces];
tag = {'East'; 'West'; 'South'; 'North'; 'Bottom'; 'Top'};
fp = fopen('G-structured.txt','w');
fprintf(fp,'cells.faces =\n');
for i=1:size(faces,1)
   fprintf(fp,' %3d %3d %3d [%s]\n', faces(i,1:3), tag{faces(i,3)});
end
nodes = [ rldecode(1:G.faces.num,diff(G.faces.nodePos), 2).' G.faces.nodes];
fprintf(fp,'\n\nfaces.nodes =\n');
fprintf(fp,' %3d %3d\n', nodes');
fprintf(fp,'\n\nfaces.neighbors =\n');
fprintf(fp,' %3d %3d\n', G.faces.neighbors');
fclose(fp);

%% UNSTRUCTURED TRIANGULAR GRID
clf;
p = sortrows([ 0.0, 1.0, 0.9, 0.1, 0.6, 0.3, 0.75; ...
               0.0, 0.0, 0.8, 0.9, 0.2, 0.6, 0.45]');
G = triangleGrid(p, delaunay(p(:,1),p(:,2)));
G = computeGeometry(G);

newplot;
plotGrid(G,'FaceColor',[0.95 0.95 0.95]); axis off;
hold on;
% centroids
text(G.cells.centroids(:,1)-0.01, G.cells.centroids(:,2)-0.01, ...
   num2str((1:G.cells.num)'),'FontSize',20);
plot(G.cells.centroids(:,1), G.cells.centroids(:,2),...
   'ok','MarkerSize',24);
% faces
text(G.faces.centroids(:,1)-0.02, G.faces.centroids(:,2)-0.01, ...
   num2str((1:G.faces.num)'),'FontSize',16);
% vertices
text(G.nodes.coords(:,1)-0.01, G.nodes.coords(:,2)-0.01, ...
   num2str((1:G.nodes.num)'),'FontSize',18);
plot(G.nodes.coords(:,1), G.nodes.coords(:,2),'sk','MarkerSize',24);
hold off;

%% Output content of cells.faces, faces.nodes, and faces.neighbors to file
fp = fopen('G-unstructured.txt','w');
faces =[ rldecode(1 : G.cells.num,diff(G.cells.facePos), 2).' G.cells.faces];
fprintf(fp,'cells.faces =\n');
fprintf(fp,' %3d %3d\n', faces');
nodes = [ rldecode(1:G.faces.num,diff(G.faces.nodePos), 2).' G.faces.nodes];
fprintf(fp,'\n\nfaces.nodes =\n');
fprintf(fp,' %3d %3d\n', nodes');
fprintf(fp,'\n\nfaces.neighbors =\n');
fprintf(fp,' %3d %3d\n', G.faces.neighbors');
fclose(fp);

%% UNSTRUCTURED VORONOI
clf;
p = [ 0.0, 1.0, 1.0, 0.0, 0.0; ...
      0.0, 0.0, 1.0, 1.0, 0.7]';
G = pebi(triangleGrid(p, delaunay(p(:,1),p(:,2))));
G = computeGeometry(G);
newplot;
plotGrid(G,'FaceColor',[0.95 0.95 0.95]); axis off;
%plotGrid(tri2grid( delaunay(p(:,1),p(:,2)), p), 'FaceColor','none','EdgeColor','r');
hold on;
% centroids
text(G.cells.centroids(:,1)-0.01, G.cells.centroids(:,2), ...
   num2str((1:G.cells.num)'),'FontSize',20);
plot(G.cells.centroids(:,1), G.cells.centroids(:,2),...
   'ok','MarkerSize',24);
% faces
text(G.faces.centroids(:,1)-0.02, G.faces.centroids(:,2)-0.01, ...
   num2str((1:G.faces.num)'),'FontSize',16);
% vertices
text(G.nodes.coords(:,1)-0.02, G.nodes.coords(:,2)-0.005, ...
   num2str((1:G.nodes.num)'),'FontSize',18);
plot(G.nodes.coords(:,1), G.nodes.coords(:,2),'sk','MarkerSize',24);
hold off;

%% Output content of cells.faces, faces.nodes, and faces.neighbors to file
fp = fopen('G-voronoi.txt','w');
faces =[ rldecode(1 : G.cells.num,diff(G.cells.facePos), 2).' G.cells.faces];
fprintf(fp,'cells.faces =\n');
fprintf(fp,' %3d %3d\n', faces');
nodes = [ rldecode(1:G.faces.num,diff(G.faces.nodePos), 2).' G.faces.nodes];
fprintf(fp,'\n\nfaces.nodes =\n');
fprintf(fp,' %3d %3d\n', nodes');
fprintf(fp,'\n\nfaces.neighbors =\n');
fprintf(fp,' %3d %3d\n', G.faces.neighbors');
fclose(fp);

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
