%% Show pillar grid
% Start with pillars, add points and mark fault pillars, add stacks of
% cells, and finally plot the grid with the fault surface marked

doprint = true;

%% Build data model
grdecl  = simpleGrdecl([4, 2, 3], 0.12, 'flat', true);
[X,Y,Z] = buildCornerPtPillars(grdecl,'Scale',true);
[x,y,z] = buildCornerPtNodes(grdecl);
G       = processGRDECL(grdecl);

%%
% Plot pillars
clf
plot3(X',Y',Z','k'); 
set(gca,'zdir','reverse'), view(35,35), axis off
axis([-0.05 1.05 -0.05 1.05 0 0.4]);
if doprint
   print -depsc2 showCPgrid-pillars.eps;
   !epstopdf showCPgrid-pillars.eps;
end

%%
% Plot points on pillars, mark pillars with faults red
hold on; I=[3 8 13]; 
hpr = plot3(X(I,:)',Y(I,:)',Z(I,:)','r','LineWidth',2);
hpt = plot3(x(:),y(:),z(:),'o');
hold off;
if  doprint
   print -depsc2 showCPgrid-pts.eps;
   !epstopdf showCPgrid-pts.eps;
end

%%
% Plot two stacks of cells
args = {'FaceColor'; 'r'; 'EdgeColor'; 'k'};
hcst = plotGrid(G,[1:8:24 7:8:24],'FaceAlpha', 0.1, args{:});
if doprint
   print -dpng showCPgrid-cells.png;
   !convert  -crop 1000x760+120+35 showCPgrid-cells.png showCPgrid-cells.png;
   !convert -resize 500x380 showCPgrid-cells.png showCPgrid-cells.eps;
end

%%
% Plot cells and fault surface
delete([hpt; hpr; hcst]);
plotGrid(G,'FaceAlpha', 0.15, args{:});
plotFaces(G,find(G.faces.tag>0),'FaceColor','b','FaceAlpha',0.4);
if doprint
   print -dpng showCPgrid-grid.png;
   !convert -crop 1000x760+120+35 showCPgrid-grid.png showCPgrid-grid.png;
   !convert -resize 500x380 showCPgrid-grid.png showCPgrid-grid.eps;
end

%% Show subdivision of fault face
clf;
% First a regular subdivision
subplot('position',[0.025 0.525 0.45 0.45]);
grdecl = simpleGrdecl([2, 1, 3], 0.12);
G  = processGRDECL(grdecl);
cellNo = rldecode(1:G.cells.num,diff(G.cells.facePos),2).';
%
cells = 1:2:6; i = false(G.cells.num,1); i(cells) = true; i = i(cellNo); 
j = false(G.faces.num, 1); j(G.faces.tag>0)=true; j = j(G.cells.faces(:,1));
plotFaces(G,G.cells.faces(i&j,1), 'FaceColor','r', 'FaceAlpha', 0.3);
%
cells = 2:2:6; i = false(G.cells.num,1); i(cells) = true; i = i(cellNo);
j = false(G.faces.num, 1); j(G.faces.tag>0)=true; j = j(G.cells.faces(:,1));
plotFaces(G,G.cells.faces(i&j,1), 'FaceColor','b', 'FaceAlpha', 0.3);
% 
plotGrid(G,'FaceColor','k','EdgeColor','k','FaceAlpha',0.025);
set(gca,'DataA',[1 1 0.5]);
set(gca,'zdir','reverse'), view(60,35), axis off, zoom(1.5);
%
subplot('position',[0.025 0.025 0.45 0.45]);
plotCellData(G,(1:6)',1:2:6,'EdgeColor','k','FaceAlpha',0.8);
set(gca,'DataA',[1 1 0.5]);
set(gca,'zdir','reverse'), view(60,35), axis off, zoom(1.5);
camdolly(0,-0.3,0);

%%
% Then an irregular subdivision
% Use same model, but perturbe some of the coordinates so that the
% intersections are no longer rectangular
subplot('position',[0.525 0.525 0.45 0.45]);
grdecl.ZCORN( 7:8:end) = grdecl.ZCORN( 7:8:end) - 0.1;
grdecl.ZCORN( 8:8:end) = grdecl.ZCORN( 8:8:end) - 0.04;
grdecl.ZCORN(11:8:end) = grdecl.ZCORN(11:8:end) - 0.05;
grdecl.ZCORN(27:8:end) = grdecl.ZCORN(27:8:end) + 0.05;
G  = processGRDECL(grdecl);
cellNo = rldecode(1:G.cells.num,diff(G.cells.facePos),2).';
%
cells = 1:2:6; i = false(G.cells.num,1); i(cells) = true; i = i(cellNo); 
j = false(G.faces.num, 1); j(G.faces.tag>0)=true; j = j(G.cells.faces(:,1));
plotFaces(G,G.cells.faces(i&j,1), 'FaceColor','r', 'FaceAlpha', 0.3);
%
cells = 2:2:6; i = false(G.cells.num,1); i(cells) = true; i = i(cellNo);
j = false(G.faces.num, 1); j(G.faces.tag>0)=true; j = j(G.cells.faces(:,1));
plotFaces(G,G.cells.faces(i&j,1), 'FaceColor','b', 'FaceAlpha', 0.3);
% 
plotGrid(G,'FaceColor','k','EdgeColor','k','FaceAlpha',0.025);
set(gca,'DataA',[1 1 0.5]);
set(gca,'zdir','reverse'), view(60,35), axis off, zoom(1.5);
%
subplot('position',[0.525 0.025 0.45 0.45]);
plotCellData(G,(1:6)',1:2:6,'EdgeColor','k','FaceAlpha',0.8);
set(gca,'DataA',[1 1 0.5]);
set(gca,'zdir','reverse'), view(60,35), axis off, zoom(1.5);
camdolly(0,-0.3,0);
%
if doprint
   print -dpng showCPgrid-faultdiv.png;
   !convert -crop 1100x760+50+0 showCPgrid-faultdiv.png showCPgrid-faultdiv.png;
   !convert -resize 550x380 showCPgrid-faultdiv.png showCPgrid-faultdiv.eps;
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
