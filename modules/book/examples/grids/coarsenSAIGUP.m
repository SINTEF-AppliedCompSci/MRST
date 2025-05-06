%% Coarsen Real Models:  SAIGUP Shallow-Marine Reservoir Model
% The <http://www.fault-analysis-group.ucd.ie/Projects/SAIGUP.html>SAIGUP
% project is a systematic assessment of uncertainty in reserves
% and production estimates within an objectively defined geological
% parameterisation encompassing the majority of European clastic oil
% reservoirs. A broad suite of shallow marine sedimentological reservoir
% types are indexed to continuously varying 3D anisotropy and heterogeneity
% levels. Structural complexity ranges from unfaulted to compartmentalised,
% and fault types from transmissive to sealing. Several geostatistical
% realisations each for the geologically diverse reservoir types covering
% the pre-defined parameter-space are up-scaled, faulted and simulated with
% an appropriate production strategy for an approximately 20 year period.
% Herein, we will inspect in detail one instance of the model, which can be
% downloaded from the <http://www.sintef.no/Projectweb/MRST>MRST webpage

%% Load dataset
mrstModule add libgeometry coarsegrid
filename = fullfile(getDatasetPath('SAIGUP'), 'SAIGUP.GRDECL');
G = processGRDECL(readGRDECL(filename));
try
   G = mcomputeGeometry(G);
catch %#ok<CTCH>
   G = computeGeometry(G);
end

%%
% We construct a coarse grid by partitioning the grid uniformly as 6x12x3
% coarse blocks in index space. This process partitions all cells in the
% logical 40x120x20 grid, including cells that are inactive.  We therefore
% compress the partitioning to remove blocks that contain no active cells,
% and then renumber the overall partitioning. Some of the blocks may
% contain disconnected cells because of erosion, faults, etc. We therefore
% postprocess the grid in physical space and split disconnected blocks.
p = partitionUI(G,[6 12 3]);
p = compressPartition(p);
p = processPartition(G,p);

%%
% We have now obtained a partitioning consisting of 243 blocks, in which
% each coarse block consists of a set of connected cells in the fine grid.
% To show the partition, we plot the coarse blocks using the explosion
% view technique
figure('Position',[0 60 800 420]);
axes('Position',[.01 .11 .775 .81]);
explosionView(G, p, 'EdgeAlpha', .1);
axis tight off; set(gca,'DataAspect',[1 1 0.1]); 
view(-65,55); zoom(1.4); camdolly(0,-0.2,0)
colormap(.3*(2*colorcube(max(p))+ones(max(p),3)));
set(gca,'Clipping','off')

%%
% Build the coarse-grid and inspect the distribution of block volumes. The
% original fine-scale model has the peculiar characteristic that all cells
% have the same size. In the coarsened model there is almost two orders
% difference between the smallest and largest blocks.
CG = generateCoarseGrid(G, p); P=p;
CG = coarsenGeometry(CG);

axes('Position',[.75 .12 .18 .8]);
h = barh(CG.cells.volumes,'b');
set(gca,'xdir','reverse','YAxisLocation','right',...
   'Fontsize',8,'XLim',[0 12]*1e6);
set(h,'FaceColor',[.7 .7 .7],'EdgeColor',[.6 .6 .6]);

%%
% To get a more even size distribution for the coarse blocks, we will
% remove some of the smallest blocks by merging them with the neighbor that
% has the smallest block volume. This is done repeatedly until the volumes
% of all blocks are above a certain lower threshold. First, we visualize all
% blocks that have volume less than ten percent of the mean block volume.
blockVols = CG.cells.volumes;
meanVol   = mean(blockVols);
show = blockVols<.1*meanVol;
figure
plotGrid(G,'FaceColor','none','EdgeAlpha',.05);
bcol = zeros(max(p),1); bcol(show)=1:sum(show); 
plotCellData(CG,bcol, show);
axis tight off; set(gca,'DataAspect',[1 1 0.1]); 
view(-65,55); zoom(1.4); camdolly(0,-0.2,0)
colormap(lines);
set(gca,'Clipping','off')

%%
% The merging algorithm is quite simple: we compute the block volumes,
% select the block with the smallest volume, and then merge this block with
% one of its neighbors, depending on some criterion (e.g., the neighbor
% with the largest or smallest block volume). Then we update the partition
% vector by relabling all cells in the block with the new block number,
% compress the partition vector to get rid of empty entries, regenerate a
% coarse grid, recompute block volumes, pick the block with the smallest
% volume in the new grid, and so on. In each iteration, we plot the
% selected block and its neighbors.

blockVols = CG.cells.volumes;
meanVol   = mean(blockVols);
[minVol, block] = min(blockVols);
while minVol<.1*meanVol
   % Find all neighbors of the block
   clist = any(CG.faces.neighbors==block,2);
   nlist = reshape(CG.faces.neighbors(clist,:),[],1);
   nlist = unique(nlist(nlist>0 & nlist~=block));

   % Alt 1: merge with neighbor having smallest volume
   % [~,merge] = min(blockVols(nlist));
   
   % Alt 2: sort neigbors by cell volume and merge with the one in the
   % middle
   % [~,ind] = sort(blockVols(nlist),1,'descend');
   % merge = ind(round(numel(ind)/2));
   
   % Alt 3: merge with neighbor having largest volume
   [~,merge] = max(blockVols(nlist));
   
   figure; 
   plotBlockAndNeighbors(CG, block, ...
      'PlotFaults', [false, true], 'Alpha', [1 .8 .8 .8]);

   % Update partition vector
   p(p==block)  = nlist(merge);
   p = compressPartition(p);
   
   % Regenerate coarse grid and pick the block with the smallest volume
   CG = generateCoarseGrid(G, p);
   CG = coarsenGeometry(CG);
   blockVols = CG.cells.volumes;
   [minVol, block] = min(blockVols);
end

%%
% Make new plot with the revised coarse grid
figure('Position',[0 60 800 420]);
axes('Position',[.01 .11 .775 .81]);
explosionView(G, CG.partition, 'EdgeAlpha', .1);
axis tight off; set(gca,'DataAspect',[1 1 0.1]); 
view(-65,55); zoom(1.4); camdolly(0,-0.2,0)
colormap(.3*(2*colorcube(CG.cells.num)+ones(CG.cells.num,3)));
set(gca,'Clipping','off')

axes('Position',[.75 .12 .18 .8]);
h = barh(CG.cells.volumes,'b');
set(gca,'xdir','reverse','YAxisLocation','right',...
   'Fontsize',8,'XLim',[0 12]*1e6);
set(h,'FaceColor',[.7 .7 .7],'EdgeColor',[.6 .6 .6]);

%%
% In the partition above, there were several blocks that were only
% connected through small face areas in their interior. To try to improve
% this problem, we recompute the the partitioning, but now disregard
% couplings across small faces when postprocessing the initial, uniform
% partition
p = partitionUI(G,[6 12 3]);
p = compressPartition(p);
p = processPartition(G, p, G.faces.areas<250);
CG1 = generateCoarseGrid(G,p);
CG1 = coarsenGeometry(CG1);

%%
% The next step is to redo the merging of small blocks. The algorithm used
% above was written to make the visualization of the changing grid as
% simple as possible. The implementation is not very efficient since we in
% each step regenerate the coarse grid and compute all the block volumes.
% This time, we will use another implementation that is a bit more
% involved, but also more efficient.

blockVols       = CG1.cells.volumes;
meanVol         = mean(blockVols);
newBlockList    = 1:max(p);
[minVol, block] = min(blockVols);
while minVol<.1*meanVol
   nlist = reshape(CG1.faces.neighbors(any(CG1.faces.neighbors==block,2),:),[],1);
   nlist = newBlockList(nlist(nlist>0));
   nlist = unique(nlist(nlist~=block));
   
   [nVol,merge] = min(blockVols(nlist));
   
   newBlockNo   = nlist(merge);
   p(p==block)  = newBlockNo;
   newBlockList(block)   = newBlockNo;
   blockVols(block)      = inf;
   blockVols(newBlockNo) = minVol + nVol;
   
   [minVol, block] = min(blockVols);
end
p = compressPartition(p);
CG1 = generateCoarseGrid(G, p);
CG1 = coarsenGeometry(CG1);

%%
% Make new plot with the revised coarse grid
figure('Position',[0 60 800 420]);
axes('Position',[.01 .11 .775 .81]);
explosionView(G, CG1.partition, 'EdgeAlpha', .1);
axis tight off; set(gca,'DataAspect',[1 1 0.1]); 
view(-65,55); zoom(1.4); camdolly(0,-0.2,0)
colormap(.3*(2*colorcube(CG1.cells.num)+ones(CG1.cells.num,3)));
set(gca,'Clipping','off')

axes('Position',[.75 .12 .18 .8]);
h = barh(CG1.cells.volumes,'b');
set(gca,'xdir','reverse','YAxisLocation','right',...
   'Fontsize',8,'XLim',[0 12]*1e6);
set(h,'FaceColor',[.7 .7 .7],'EdgeColor',[.6 .6 .6]);

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
