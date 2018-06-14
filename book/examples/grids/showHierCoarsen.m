%% Hierarchical Coarsening
% This is a relatively simple example that illustrates the basic idea
% behind hierarchical coarsening. The model has two geological properties,
% unit and lithofacies assemblage (LFA), that are used to generate
% permeability. These two partitions are applied recursively together with
% a nested set of uniform partitions with permeability as an indicator.

mrstModule add coarsegrid agglom

% Load facies data
exdir = fullfile(mrstPath('agglom'), 'examples');
imload = @(fn) ...
   flipud(double(sum(imread(fullfile(exdir, 'data', fn)), 3))) .';
f  = imload('facies1.png'); f(2,16)=max(f(:));
pl = 3-compressPartition(f(:) + 1);


% Create grid
G  = computeGeometry(cartGrid(size(f))); clear f;
x = G.cells.centroids;

% Make unit data
pu = (x(:,2)>15+0.25*x(:,1))+1;

% Make fault blocks
faults = [12:42:29*42 1692:41:2840 432:41:1630];
pf = processPartition(G, ones(G.cells.num,1), faults);

% Generate and plot a simple permeability model
rng(1000);
K  = 1 + (pl-1)*8+.3*pu.*randn(G.cells.num,1); I = K-min(K)+.1;
figure(1); clf, 
plotCellData(G,I,'EdgeColor','none');
axis equal off; colormap(.8*jet(128)+.2*ones(128,3))


%% Apply the geological features recursively
% In addition to the geological features, we create a hierarchy of refined
% partitions that will be used to further subdivide the high-permeable LFA
pc = ones(G.cells.num,3);
for i=1:4
    pc(:,i) = partitionCartGrid(G.cartDims,[5 5].*2^(i-1));
end
p = applySuccessivePart(processPartition(G,pu,faults),G, I, 8, [pl pc]);

% Visualize results
figure(2); clf, set(gcf,'Position',[670 490 760 470]);
subplot(2,3,1)
plotCellData(G, pl, 'EdgeColor','none'); axis equal off
subplot(2,3,2);
plotCellData(G, pu, 'EdgeColor','none'); axis equal off
subplot(2,3,3)
plotCellData(G, pl+2*(pu-1), 'EdgeColor','none'); axis equal off
plotFaces(G, faults,'EdgeColor','k','LineWidth',1)
subplot(2,3,5); cla
plotCellData(G, pl+2*(pu-1),'EdgeColor','none'); axis equal off
outlineCoarseGrid(G, p, 'k','LineWidth',1);
colormap([.5 .5 1; .2 .2 .8; .3 .6 .6; .1 .4 .4]); 

%% Merge small blocks outside of high-permeable lithofacies
% We merge all blocks that have three cells or less in the first LFA, which
% is assumed to have less permeability than the second lithofacies. To this
% end, we remove cell connections within LFA two as well as
% connections across different units, lithofacies, and faults.
subplot(2,3,6), cla
plotCellData(G, pl+2*(pu-1),'EdgeColor','none'); axis equal off

qb = [0; processPartition(G,pl+2*(pu-1), faults)];
ql = [0; pl];
T = 2*(qb(G.faces.neighbors(:,1)+1)==qb(G.faces.neighbors(:,2)+1))-1;
T(ql(G.faces.neighbors(:,1)+1)==2) = -1;
pm = mergeBlocksByConnections(G, p, T, 4);

outlineCoarseGrid(G,pm, 'k','LineWidth',1);

for i=[1:3 5:6], 
    subplot(2,3,i), 
    set(gca,'Position',get(gca,'Position')+[-.025 -.025 .05 .05]);
end
