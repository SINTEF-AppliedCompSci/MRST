%% Examples of How to Generate Coarse Grids
% We show three examples: a 2x2 partition of a 4x4 Cartesian grid, a 2x2
% partition of a facies model with subpartition of coarse faces, and a 3D
% model with subdivision of coarse faces

mrstModule add coarsegrid;

%% First example: 2x2 partition of a 4x4 Cartesian grid
G = computeGeometry(cartGrid([4,4]));
p = partitionUI(G, [2,2]);
CG = generateCoarseGrid(G, p);

clf
plotCellData(G,(1:G.cells.num)','LineStyle',':');
plotGrid(CG,'faceColor','none','LineWidth',3);
colormap((jet(G.cells.num)+repmat(2,G.cells.num,3))/3);

%% 
% Show cell/block indices
CG = coarsenGeometry(CG);
tg = text(G.cells.centroids(:,1), G.cells.centroids(:,2), ...
   num2str((1:G.cells.num)'),'FontSize',20, 'HorizontalAlignment','center');
tcg = text(CG.cells.centroids(:,1), CG.cells.centroids(:,2), ...
   num2str((1:CG.cells.num)'),'FontSize',24, 'HorizontalAlignment','center');
axis off;
set(tcg,'BackgroundColor','w','EdgeColor','none');

%%
% Show indices of the faces in the fine/coarse grids
delete([tg; tcg]);
tg = text(G.faces.centroids(:,1), G.faces.centroids(:,2), ...
   num2str((1:G.faces.num)'),'FontSize',18, 'HorizontalAlignment','center');
tcg = text(CG.faces.centroids(:,1), CG.faces.centroids(:,2), ...
   num2str((1:CG.faces.num)'),'FontSize',24, 'HorizontalAlignment','center');
set(tcg,'BackgroundColor','w','EdgeColor','k');

%% Second example: 2D face partition 
clf;
G  = computeGeometry(cartGrid([8, 8], [1 1]));
f  = @(c) sin(3*pi*(c(:,1)-c(:,2)));
pf = 1 + (f(G.cells.centroids) > 0);

plotCellData(G, pf,'edgecolor','none'); axis off;
colormap(.2*gray(2)+.8*ones(2,3));

pv = partitionCartGrid(G.cartDims, [2 2]);
pf = cellPartitionToFacePartition(G,pf);
pf = processFacePartition(G, pv, pf);
CG = generateCoarseGrid(G, pv, pf);

CG = coarsenGeometry(CG);
cmap = lines(CG.faces.num);
for i=1:CG.faces.num
    plotFaces(CG,i,'LineWidth',6,'EdgeColor', cmap(i,:));
end
text(CG.faces.centroids(:,1), CG.faces.centroids(:,2), ...
        num2str((1:CG.faces.num)'),'FontSize',20,'HorizontalAlignment','center');


%% Third example: 3D face partition
clf
G = computeGeometry(cartGrid([20 20 6]));
c = G.cells.centroids;
G = removeCells(G, ...
   (c(:,1)<10) & (c(:,2)<10) & (c(:,3)<3));
plotGrid(G,'FaceColor',[1 1 .7]); view(3); axis off

%%
clf
p = partitionUI(G,[2, 2, 2]);
q = partitionUI(G,[4, 4, 2]);
CG = generateCoarseGrid(G, p, ...
   cellPartitionToFacePartition(G,q));

plotCellData(CG,(1:max(p))');
plotFaces(CG,1:CG.faces.num,...
   'FaceColor' , 'none' , 'LineWidth' ,2);
view(3); axis off
colormap(.5*(jet(128)+ones(128,3)));

%%
% Add face centroids
CG = coarsenGeometry(CG);
plotCentroids = @(pts, varargin) ...
   plot3(pts(:,1), pts(:,2), pts(:,3), varargin{:});
hold on
h=plotCentroids(CG.faces.centroids, 'k*');
hold off;

%%
% The coarse grid also contains lookup tables for mapping blocks and
% interfaces in the coarse grid to cells and faces in the fine grid. To
% demonstrate this, we visualize a single coarse face consisting of several
% fine faces along with the cells its block neighbors in red and blue
% respectively on the fine grid.
face  = 66;
sub   = CG.faces.connPos(face):CG.faces.connPos(face+1)-1;
ff    = CG.faces.fconn(sub);
neigh = CG.faces.neighbors(face,:);

clf
show = false(1,CG.faces.num);
show(boundaryFaces(CG)) = true;
show(boundaryFaces(CG,neigh)) = false;
plotFaces(CG, show,'FaceColor',[1 1 .7]);
plotFaces(CG,boundaryFaces(CG,neigh),'FaceColor','none','LineWidth', 2);
plotFaces(G, ff, 'FaceColor', 'g')
plotGrid(G, p == neigh(1), 'FaceColor', 'none', 'EdgeColor', 'r')
plotGrid(G, p == neigh(2), 'FaceColor', 'none', 'EdgeColor', 'b')
view(3); axis off

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
