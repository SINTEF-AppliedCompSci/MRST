%% More About Coarsegrid Geometry/Topology
% In this example, we will take a closer look at ways of working with the
% geometry and topology of the coarse grids

mrstModule add coarsegrid
%% Subdivision of Coarse Faces
% So far, we have assumed that there is only a single connection between
% two neighboring coarse blocks. This connection is built up of all faces
% between pairs of fine cells that belong to the two different coarse
% blocks. However, we can make coarse grids that have multiple connections
% between neighboring grid blocks, as the following examples shows.
G = computeGeometry(cartGrid([20 20 6]));
c = G.cells.centroids;
G = removeCells(G, (c(:,1)<10) & (c(:,2)<10) & (c(:,3)<3));

p = partitionUI(G,[2, 2, 2]); p = compressPartition(p); mp = max(p);
q = partitionUI(G,[4, 4, 2]);
CG = generateCoarseGrid(G, compressPartition(p), ...
                        cellPartitionToFacePartition(G,q));
clf
plotCellData(CG,(1:mp)');
plotFaces(CG,1:CG.faces.num,'FaceColor' , 'none' , 'LineWidth' ,2);
view(3); axis off
colormap(.5*(jet(128)+ones(128,3)));


%% Adding geometry information to the coarse grid
% Unlike the fine grids, the coarse grid structure only contains
% topological information by default. This means, in particular, that
% coarse grids do not have nodes. By calling coarsenGeometry, we can get
% coarse block and face centroids, volumes, and so on.
CG = coarsenGeometry(CG);

plotPts = @(pts, varargin) plot3(pts(:,1), pts(:,2), pts(:,3), varargin{:});
hold on
h=plotPts(CG.faces.centroids, 'k*');
hold off;


%%
% Let us look more in detail on the relationship between centroids in the
% fine and coarse grid
i = false(mp,1); i(4) = true;
cla, hold
cg_cent = CG.cells.centroids(i,:);   plotPts(cg_cent, 'ok','MarkerFaceColor',[.5 .5 .5]);
g_cent = G.cells.centroids(p==4,:);  plotPts(g_cent, '.');
plotCellData(CG,(1:mp)',~i);
plotFaces(CG,gridCellFaces(CG,find(~i)),'FaceColor','none','LineWidth',2);
view(-35,15);
legend({'Coarse centroids', 'Fine Centroids'}, 'Location', 'SouthOutside')

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

% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
