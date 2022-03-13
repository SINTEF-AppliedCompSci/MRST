%% More about partition vectors
% In this tutorial, we will take a closer look at partition vectors and
% discuss how different types of partitions can be combined into one.
mrstModule add coarsegrid

%% Partition according to polar coordinate
% In the first example, we partition the box [-1,1]x[-1,1] into nine
% different blocks using the polar coordinates of the cell centroids. The
% first block is defined as r<=.3, while the remaining eight are defined by
% segmenting 4 theta/pi:
clf
G = cartGrid([11, 11],[2,2]);
G.nodes.coords = bsxfun(@minus, G.nodes.coords, 1);
G = computeGeometry(G);
c = G.cells.centroids;
[th,r] = cart2pol(c(:,1),c(:,2));
p = mod(round(4*th/pi)+4,4)+1;
p(r<.3) = max(p)+1;

plotCellData(G,p,'EdgeColor',[.7 .7 .7]);
outlineCoarseGrid(G,p,'k');
caxis([.5 max(p)+.5]);
colormap(.5*(jet(max(p))+ones(max(p),3)));
set(colorbar,'YTick',1:max(p),'FontSize',16);
axis off

%%
% In the figure it is simple to distinguish nine different coarse blocks,
% but since we only have five colors, there are actually five blocks. Four
% of these are disconnected in the sense that they contain cells that
% cannot be connected by a path that only crosses faces between cells
% inside the block.  We use the routine 'processPartition' to split these
% blocks
p = processPartition(G, p);
figure
plotCellData(G,p,'EdgeColor',[.7 .7 .7]);
outlineCoarseGrid(G,p,'k');
caxis([.5 max(p)+.5]);
colormap(.5*(jet(max(p))+ones(max(p),3)));
set(colorbar,'YTick',1:max(p),'FontSize',16);
axis off


%% No need for logical indices
% As you have seen in the previous example, the partition vector is just a
% vector that tells the index of the block a given cell belongs to. This
% means that we can generate an arbitrary coarse grid without using any
% information relating to the geometry or topology of the grid.
%
% For instance, if we want to partition a grid uniformly in space based on
% cell coordinates, the product of two periodic functions will do nicely.
% To illustrate this, we can divide a rectangular grid into a 3x3 coarse
% grid by exploiting that the sine function changes sign in intervals of
% pi.
nx = 3; ny = 3;
G = computeGeometry(cartGrid([20, 20], [1 1]));

% First create the periodic function
f = @(v) sin(pi*nx*v(:,1)) .* sin(pi*ny*v(:,2));
% Evaluate it in the centroids
fval = f( G.cells.centroids);

% We divide the grid into two parts based on the sign of the function
% and postprocess to create connected domains.
p = double(fval > 0) + 1;
p = processPartition(G, p);

% This vector can generate coarse grids just as partitionUI did.
CG = generateCoarseGrid(G, p);
clf
subplot(1,2,1); title('f(x)')
v = reshape(f(G.cells.centroids), G.cartDims(1), G.cartDims(2));
surf(v)

subplot(1,2,2); title('Resulting partition')
plotCellData(G, mod(p, 13), 'EdgeColor', 'w');
colormap(jet)

%% Combining multiple partition vectors
% A very efficient way to generate coarse partitions is to combine two or
% more partition vectors that each implement a coarsening principle. As an
% examples, let us consider a case where our geological model contains two
% different facies. We can use this information, together with a standard
% partition, to make a coarse grid that adapts to the geology in the sense
% that we only get coarse blocks that contain a single facies. This could,
% for instance, be important when upscaling properties that are very
% different.
clf

% Make a diagonal and periodic facies distribution
f  = @(c) sin(4*pi*(c(:,1)-c(:,2)));
pf = 1 + (f(G.cells.centroids) > 0);
subplot(1,2,1);plotCellData(G,pf); view(2);

% Make a rectangular coarse partition
pc = partitionCartGrid(G.cartDims, [4 4]);
subplot(1,2,2), plotCellData(G,pc); view(2); colormap(jet(128));

%%
% Alternative 1: Collect the partitions as columns in an array and the use
% 'unique' to find all unique combinations. These will give the
% intersection of the two. This technique can be generalized to an
% arbitrary number of partition vectors, but is costly since 'unique'
% performs a search.
tic, [b,~,p] = unique([pf, pc], 'rows' ); toc

%%
% Alternative 2: We can treat the partition vectors as multiple subscripts,
% and use this to compute a linear index. This index will not be contiguous
% and thus needs to be compressed. This method is fast, but more difficult
% to extend to multiple partition vectors using a single statement.
tic, q = compressPartition(pf + max(pf)*pc); toc

%%
% Plot the resulting grid
figure
plotCellData(G,p);
outlineCoarseGrid(G, p, 'k','LineWidth',2);
axis off
mp = max(p);
colormap(.5*(colorcube(mp)+ones(mp,3)))

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
