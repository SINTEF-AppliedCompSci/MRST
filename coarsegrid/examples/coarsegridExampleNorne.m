%% Coarsening the Norne Simulation Model
% We show how to coarsen the simulation model of the Norne field by
% partitioning it uniformly in logical Cartesian space. We visualize some
% of the coarse blocks and show how they are connected with their
% neighbors.

mrstModule add coarsegrid

%% Read and process model
% The Norne model is one of the standard data sets that have been used a
% lot with MRST. Unfortunately, you have to download the model manually.
% How to do this, can be found by typing |getDatasetInfo('norne')|. We pick
% the part of the model that represents the main reservoir and disregard
% the stack of twelve cells that are disconnected from the rest of the
% model.
dpath  = getDatasetPath('norne');
grdecl = fullfile(dpath, 'NORNE.GRDECL');
grdecl = readGRDECL(grdecl);
G = processGRDECL(grdecl);
G = computeGeometry(G(1));

clf,
plotGrid(G,'EdgeAlpha',.1,'FaceColor',[.9 .9 .7]); 
view(70,70), zoom(2), axis tight off
set(gca,'Zdir','normal'); camlight headlight; set(gca,'Zdir','reverse');

%% Partition the grid in logical space
% To construct a coarse grid, we try to partition the grid uniformly as
% 6x12x3 coarse blocks in index space. This will partition all cells in the
% logical 46x112x22 grid, including cells that are inactive. The number of
% active cells within each coarse block is shown in the bar plot below.

p = partitionUI(G,[6 12 3]);
newplot, subplot(3,1,1)
bar(accumarray(p,1)); set(gca,'XLim',[0 225]);
title('Unprocessed');

%%
% As we can see from the bar plot, there are several coarse block that
% contain no active cells. We therefore postprocess the partitioning to
% remove blocks that contain no active cells, and then renumber the overall
% partitioning, giving a new total of 168 blocks.
p = compressPartition(p); max(p)
subplot(3,1,2)
bar(accumarray(p,1)); set(gca,'XLim',[0 225]);
title('Compressed');

%%
% Because the partitioning has been performed in logical index space, we
% have so far disregarded the fact the some of the blocks may contain
% disconnected cells because of erosion, faults, etc. We therefore
% postprocess the grid in physical space and split disconnected blocks.
p = processPartition(G,p); m=max(p);
subplot(3,1,3)
bar(accumarray(p,1)); set(gca,'XLim',[0 225]);
title('Processed');

assert (all(accumarray(p, 1) > 0))

%%
% We have now obtained a partitioning consisting of 223 blocks, in which
% each coarse block consists of a set of connected cells in the fine grid.
% To show the partitioning, we plot the coarse blocks using a random and
% cyclic color scheme for the blocks.
newplot
subplot('position',[0.025 0.025 0.95 0.95])
plotCellData(G,p,'EdgeAlpha',.5);
axis tight off; view(0,75), colormap(colorcube(m));

%%
% From the plot above, it is not easy to see the shape of the individual
% coarse blocks. In the next section, we will therefore show some examples
% of how individual blocks can be visualized.


%% Build the coarse-grid
% Having obtained a partition we are satisfied with, we build the
% coarse-grid structure. This structure consists of three parts:
%
% * the cell structure giving the number of blocks and the indices of the
% cells contained in each block
% * the face structure giving the number of coarse faces and the indices of
% the neighbouring blocks
% * a cellFaces array as in the fine-grid structure
CG = generateCoarseGrid(G, p);
CG             %#ok<NOPTS> (intentional display)
CG.cells       %  (intentional display)
CG.faces       %  (intentional display)

%%
% Let us now use CG to inspect some of the blocks in the coarse grid. To
% this end, we arbitrarily pick a few blocks and inspect these block and
% their neighbours. For the first block, we plot the cells and the faces
% that have been marked as lying on a fault
clf, plotBlockAndNeighbors(CG, 55), view(-90,70)

%%
% For the second block, we only plot the cells and not the faulted faces
clf
plotBlockAndNeighbors(CG, 19, 'PlotFaults', false([2, 1]))
view(90, 70)

%%
% The third set of neighboring blocks contains more faults
clf, plotBlockAndNeighbors(CG, 191), view(100, 30)

%%
% We end the example by highlight six representative blocks, including the
% three blocks we inspected above. Notice that this way of visualization
% only uses the fine grid and the partition vector, and thus does not
% require that the coarse-grid structure has been built.
clf
blocks = [3, 19, 191, 34, 196, 55];
col = ['b', 'g', 'r', 'c', 'm', 'y'];
axes('position', [0.01, 0.25, 0.99, 0.75]);
plotGrid(G, 'EdgeColor', [0.75, 0.75, 0.75], 'FaceColor', 'w');
outlineCoarseGrid(G, p, 'FaceColor', 'none', 'LineWidth', 2);

for i = 1 : 6,
   plotGrid(G, p == blocks(i), 'FaceColor', col(i));
end
axis tight off, view(10, 90)

% Plot the chosen 6 coarse blocks
for i = 1 : 6,
   axes('position', [(i-1)/6, 0.02, 1/6, 0.25]);

   plotGrid(G, p == blocks(i), 'FaceColor', col(i));

   axis tight off, view(0,75), zoom(1.2)
end

%%
displayEndOfDemoMessage(mfilename)

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
