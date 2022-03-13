%% Second Coarse-Grid Tutorial: Partitioning More Complex Grids
% In this tutorial, we will take a look at how to partition grids that have
% more complex geometries than the simple 2D grids we used in the first
% tutorial.

mrstModule add coarsegrid;

%% Partition a Non-rectangular 2D Grid
% If the grid indices do not form a logical rectangle, we need to use the
% function 'partitionUI', which computes a uniform partition of the
% bounding box of the grid. To illustrate, we consider a rectilinear grid
% with a semi-circular cut-out along the lower edge.
dx = 1-0.5*cos((-1:0.1:1)*pi); 
x = -1.15+0.1*cumsum(dx);
y = 0:0.05:1;
G = computeGeometry(tensorGrid(x, y.^.75));
r = sqrt(sum(G.cells.centroids.^2,2));
G = extractSubgrid(G, r>0.6);

clf
CG = generateCoarseGrid(G, partitionUI(G,[5 4]));
plotCellData(CG,(1:CG.cells.num)','EdgeColor','w','EdgeAlpha',.2);
plotFaces(CG,(1:CG.faces.num)', 'FaceColor','none','LineWidth',2);
colormap(.5*(colorcube(20) + ones(20,3))); axis off

%% Partitioning a Cup-Formed Grid
% 3D Cartesian grid can be partitioned in the exact same way as 2D
% rectangular grids, using partitionCartGrid for grids that completely fill
% a hexahedron and partitionUI for grids that do not. As an illustration,
% we will partition a domain formed as a cup
clf
x = linspace(-2,2,41);
G = tensorGrid(x,x,x);
G = computeGeometry(G);
c = G.cells.centroids;
r = c(:,1).^2 + c(:,2).^2+c(:,3).^2;
G = removeCells(G, (r>1) | (r<0.25) | (c(:,3)<0));
plotGrid(G); view(15,60); axis tight off

%%
% Since the grid only covers subregion of its bounding box, the partition
% vector will not be contiguous. To avoid having to treat special cases
% arising from non-contigous partition vectors, most routines in MRST
% require partition vectors to be contiguous. To get a contigous partition
% vector, we use the function 'compressPartition'
clf
p = partitionUI(G,[5 5 4]);
subplot(2,1,1); bar(accumarray(p,1)); shading flat
q = compressPartition(p);
subplot(2,1,2); bar(accumarray(q,1)); shading flat
set(gca,'XLim',[0 100]);

%%
% And then we are ready to make the coarse grid. The figure below shows two
% different ways of visualizing the grid
CG = generateCoarseGrid(G, q);

subplot(1,2,1)
explosionView(G,q,.4);
view(15,60); axis tight off
colormap(colorcube(max(q)));

subplot(1,2,2)
plotCellData(CG,(1:CG.cells.num)','EdgeColor','w','EdgeAlpha',.2);
plotFaces(CG,(1:CG.faces.num)', 'FaceColor','none','LineWidth',2);
view(15,60); axis tight off

%% Partitioning an Unstructured Grid
% The same principle applies to unstructured grids. So let us try to
% partition a Voronoi grid into a semi-rectangular coarse grid. The
% function partitionUI only works for fine grids with an underlying
% Cartesian topography, so we can rely on the function 'sampleFromBox'
% which uses the cell centroids to sample from a uniform Cartesian grid
% covering the bounding box of the grid
[x,y] = meshgrid(linspace(0,1,25));
pt = twister([x(:) y(:)]);
G = pebi( triangleGrid(pt, delaunay(pt)));
G = computeGeometry(G);
r = sqrt(sum(G.cells.centroids.^2,2));
G = extractSubgrid(G, r>0.6);
p = sampleFromBox(G, reshape(1:16,4,4));
clf
plotCellData(G,p); colormap(colorcube(16)); axis tight off

%% Create a faulted grid and partition it
% We create a grid and partition it logically in ij-space and along
% specific layers along the k-space.
grdecl = simpleGrdecl([10 10 7], .15);
G = computeGeometry( processGRDECL(grdecl) );
% Layer 1 will go from 1 to 3 - 1, layer 2 from 3 to 6 - 1 and so on
L = [1 3 6 8];
% We can easily find the thickness of the layers
diff(L)

% The partition is disconnected across the fault. processPartition can
% amend this by adding new coarse blocks.
p_initial = partitionLayers(G, [5,5], L);
p = processPartition(G, p_initial); mp = max(p);

figure
CGi = generateCoarseGrid(G,p_initial);
plotCellData(G, p_initial, 'EdgeAlpha', .2);
plotFaces(CGi,(1:CGi.faces.num),'FaceColor','none','EdgeColor','w','LineWidth',2);
title(['Before processPartition (Total blocks: ' num2str(max(p_initial)) ')'])
view(60,65); colormap(colorcube(mp)); caxis([1 mp]); axis tight off

figure
CG = coarsenGeometry( generateCoarseGrid(G,p) );
plotCellData(G, p, 'EdgeAlpha', .2);
plotFaces(CG,(1:CG.faces.num),'FaceColor','none','EdgeColor','w','LineWidth',2);
title(['After processPartition (Total blocks: ' num2str(mp) ')'])
view(60,65); colormap(colorcube(mp)); caxis([1 mp]); axis tight off

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
