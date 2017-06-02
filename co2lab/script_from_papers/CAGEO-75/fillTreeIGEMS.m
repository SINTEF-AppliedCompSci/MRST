%% Filling of a trapping tree
% In this example, we will use percolation type analysis to identify trees
% and then gradually fill the largest tree to visualize the movement of CO2
% injected at an infinitesimal rate.
mrstModule add co2lab;

%% Download data if necessary
idir = fullfile(mrstPath('co2lab'), 'data', 'igems');
if ~exist(fullfile(idir,'one_of_each'),'dir')
   disp(' -> Download data from: http://www.nr.no/IGEMS')
   disp(['    Putting data in ', idir]);
   unzip('https://www.nr.no/sites/default/files/files/one_of_each.zip', idir);
end

%% Get the grid
idir = fullfile(idir,'one_of_each');
Gt = topSurfaceGrid(...
  readIGEMSIRAP('1_OSSNP1.irap', [], 'dir', idir, 'save', false));

%% Find trap volumes, traps and trees
res = trapAnalysis(Gt, true);
trapvols = volumesOfTraps(Gt, res, []);
trees = maximizeTrapping(Gt, 'res', res, 'removeOverlap', true);

%% Plot the surface and the spill paths
% Find the longests tree
treelengths = arrayfun(@(x) numel(x.traps), trees);
[tmp, i]    = max(treelengths);
tree        = trees(i);

% Compute gradient of surface elevation and show it using a nonlinear
% graymap to indicate the relief in the top surface
[zx,zy] = gradient(reshape(Gt.cells.z,Gt.cartDims));
zgrad   = sqrt(zx.^2 + zy.^2);
plotCellData(Gt, zgrad(:),'EdgeColor','none');
colormap(flipud(gray(128)).^4);
axis equal tight off, view(-90,90);

% Plot spillpaths
cell_lines = res.cell_lines;
cc = [Gt.cells.centroids, Gt.cells.z];
[x,y,z] = deal([]);
for i=1:numel(cell_lines);
   if ~ismember(i, tree.traps) || isempty(cell_lines{i})
        continue
    end
    cl = unique(cell_lines{i}{1});
    cl = cl(res.traps(cl) == 0);
    x = [x; cc(double(cl),1); NaN]; %#ok
    y = [y; cc(double(cl),2); NaN]; %#ok
    z = [z; cc(double(cl),3); NaN]; %#ok
end
hold on
plot3(x,y,z+10, 'r', 'linewidth', 2)
hold off


%% Gradually fill the traps in the largest tree
ht = [];
for v=0:0.25:1
   injectedVolume = tree.value*v;
   fill = findPercolationVolumes(tree.traps, injectedVolume, trapvols);

   delete(ht); ht = [];
   for i=1:numel(tree.traps)
      subset = res.traps==tree.traps(i);
      h = plotGrid(Gt, subset, 'EdgeColor','none',...
         'FaceColor',[1-fill(i) 1-0.25*fill(i) .4]);
      ht = [ht; h];
   end
   drawnow; title(['Fill factor: ' num2str(v)]); pause(1);
   % print('-r300', '-dpng', ['igems-treefill-' num2str(v/.25) '.png']);
end

