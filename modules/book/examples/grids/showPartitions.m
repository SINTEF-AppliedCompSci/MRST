%% How to Partition Grids
% In this script, we go through several examples that demonstrate various
% methods for partitioning grids using MRST
mrstModule add coarsegrid

%% Partition a Cartesian 2D grid
% We use partitionUI which exploits the logical structure and creates a
% uniform grid in logical space.

figure
G = cartGrid([7,7]);
p = partitionUI(G, [2,2]);
plotCellData(G, p, 'EdgeColor', 'y');
outlineCoarseGrid(G, p, 'k');
axis tight off, 
caxis([.5 max(p)+.5]); 
colormap(0.5*(lines(max(p))+ones(max(p),3)));
set(colorbar,'YTick',1:max(p));

%% Partition a 3D grid in much the same manner
figure
G = cartGrid([10,10,4]);
p = partitionUI(G, [3,3,2]);

plotCellData(G, p, 'Edgecolor', 'w');
outlineCoarseGrid(G, p, 'EdgeColor','k','lineWidth',4);
colormap(.5*(colorcube(max(p)) + ones(max(p),3)));
view(3);
axis off

%% Partition according to polar coordinate
figure
G = cartGrid([11, 11],[2,2]);
G.nodes.coords = ...
   bsxfun(@minus, G.nodes.coords, 1);
G = computeGeometry(G);
c = G.cells.centroids;
[th,r] = cart2pol(c(:,1),c(:,2));
p = mod(round(th/pi*4)+4,4)+1;
p(r<.3) = max(p)+1;

% Plot partition
plotCellData(G,p,'EdgeColor',[.7 .7 .7]);
outlineCoarseGrid(G,p,'k');
caxis([.5 max(p)+.5]);
colormap(.5*(jet(max(p))+ones(max(p),3)));
set(colorbar,'YTick',1:max(p),'FontSize',16);
axis off

% Split blocks that are multiply connected into a set of singly connected
% blocks and plot the new partition
p = processPartition(G, p);
figure
plotCellData(G,p,'EdgeColor',[.7 .7 .7]);
outlineCoarseGrid(G,p,'k');
caxis([.5 max(p)+.5]);
colormap(.5*(jet(max(p))+ones(max(p),3)));
set(colorbar,'YTick',1:max(p),'FontSize',16);
axis off


%% Combine a facies and a Cartesian partition
clear
G = cartGrid([20, 20], [1 1]);
G = computeGeometry(G);

f  = @(c) sin(4*pi*(c(:,1)-c(:,2)));
pf = 1 + (f(G.cells.centroids) > 0);
pc = partitionCartGrid(G.cartDims, [4 4]);

% Alternative 1:
[b,~,p] = unique([pf, pc], 'rows' );

% Alternative 2:
q = compressPartition(pf + max(pf)*pc);                         %#ok<NASGU>

% Plot results
figure
plotCellData(G,p);
outlineCoarseGrid(G, p, 'k','LineWidth',2);
axis off
colormap(.5*(jet(64)+ones(64,3)))

%% Partition a cup-formed grid
figure
x = linspace(-2,2,41);
G = tensorGrid(x,x,x);
G = computeGeometry(G);
c = G.cells.centroids;
r = c(:,1).^2 + c(:,2).^2+c(:,3).^2;
G = removeCells(G, (r>1) | (r<0.25) | (c(:,3)<0));
plotGrid(G); view(15,60); axis tight off

% Make the partition vector contiguous:
figure
p = partitionUI(G,[5 5 4]);
subplot(2,1,1); bar(accumarray(p,1)); shading flat
q = compressPartition(p);
subplot(2,1,2); bar(accumarray(q,1)); shading flat
set(gca,'XLim',[0 100]);

% Visualize the partition using an exploding view
figure
explosionView(G,q,.4);
view(15,60); axis tight off
colormap(colorcube(max(q)));

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
