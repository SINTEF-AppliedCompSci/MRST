mrstModule add coarsegrid mrst-gui

%% Extract a subset of the Utsira formation for analysis
coarsefactor = 1;
[grdecl, d, petrodata] = getAtlasGrid('utsirafm', 'coarsening', coarsefactor);
pd = petrodata{1};

G = processGRDECL(grdecl{1});
G = computeGeometry(G(1));
cdims = G.cartDims;
% G = extractSubgrid(G, find(G.cells.centroids(:, 2) > 6.62e6));
% G.cartDims = cdims;
% G = computeGeometry(G);

Gt = topSurfaceGrid(G);
rock = struct('poro', repmat(pd.avgporo, G.cells.num, 1), ...
              'perm', repmat(pd.avgperm, G.cells.num, 1));
%% Compute trapping analysis and find the sorted injection trees
% An injection tree corresponds to a root node and all traps that are
% upstream from that node. By default, maximize trapping finds all
% leafnodes while accounting for overlap between the trees.
%
% We find both the trap tree with and without overlap. Once plotted, it is
% obvious that the calculated volumes with overlap are quite different from
% those without: As several leaf nodes may be downstream from the larger
% traps in any given tree, many trees will largely consist of the same
% traps. When optimizing injection with a single injection point the
% distinction is not important, but the overlap must be considered when

res = trapAnalysis(Gt, 'cell');
[trees, trapvols] = maximizeTrapping(Gt, 'res', res, 'removeOverlap', true);
trees_nooverlap = maximizeTrapping(Gt, 'res', res, 'removeOverlap', false);
%%
% several injectors are being considered.
close all
fp = get(0, 'DefaultFigurePosition');
h = figure('position', fp.*[1 1 1 .5]);

[tv, tv_noop] = deal(nan(max(res.traps), 1));

% Index into total trap index based on leaf node
tv(vertcat(trees.root)) = vertcat(trees.value);
tv_noop(vertcat(trees_nooverlap.root)) = vertcat(trees_nooverlap.value);

% Sort by the overlapping values
[tv, ind] = sort(tv, 'descend');
tv_noop = tv_noop(ind);

% Strip nan values
tv_noop(isnan(tv)) = [];
tv(isnan(tv)) = [];

bar([tv, tv_noop- tv], 'stacked')
set(gcf, 'Renderer', 'painters')
ylabel('Tree volume (m^3)')
xlabel('Tree index')
caxis( [.5 2])
set(gca, 'Yscale', 'log')
axis auto
ylim([min(tv), 1.1*max(tv)])
legend({'Without overlap', 'Overlap'},'location', 'southwest')

%% Find the ideal injection spot for each trap region
% Use the largest z value to find the ideal injection spot: This
% corresponds to the deepest possible place in the reservoir surface to
% inject.
%
% We show the ten best root notes, which maximizes stored CO2 where five
% injectors are to be placed. We also find the best single injection cell
% if we assume that the boundary between two different trap regions spills
% over to both the neighboring trap regions. 
if ishandle(h); close(h); end;
h = figure('position', fp.*[1 1 2.5 1.5]);

N = 10;
% Work out the lowest point in each root trap
trapcells = zeros(N, 1);
for i = 1:N
    r = res.trap_regions == trees(i).root;
    region = find(r);
    tmp = Gt.cells.z;
    tmp(~r) = -inf;
    
    [minpoint, cell] = max(tmp);
    trapcells(i) = cell;
end
% Get the best cell on a ridge
[bestSingleCell, largestVol, allFaces, point, traps] = findOptimalInjectionPoint(Gt, res);

% Vast differences in values, use log10
largest = log10(largestVol);
v       = log10(tv);
treestart = [largest; v(1:N)];
trapcells = [bestSingleCell; trapcells];

clf;
% Plot all traps colorized by their largest tree volume
plotted = false(G.cells.num, 1);
for i = 1:numel(trees)
    cells = ismember(res.traps, trees(i).traps);
    plotCellData(G, repmat(v(i), G.cells.num, 1), cells & ~plotted);
    plotted = plotted | cells;
end

plotGrid(G, res.traps == 0, 'edgec', 'none', 'facec', [1 1 1]*.9)

plotBar = @(data, cells, format) plotGridBarchart(G, data, cells, ...
                            'widthscale', 5*(1/coarsefactor), 'facecolor', [], ...
                            'fontsize', 15, ...
                            'EdgeColor', 'none', 'textformat', format);
% Setup function handle to show percent of total trapped volume
perc = @(x) 100*(10^x)/sum(tv);
plotBar(treestart, trapcells, @(x) sprintf('%1.1f%%', perc(x)))

outlineCoarseGrid(G, res.trap_regions, 'facecolor', 'none')
% set(gca,'ZDir','normal');
% camlight head
% lighting phong
% set(gca,'ZDir','reverse');

axis tight off
view(64, 72)
caxis([min(treestart), max(treestart)]);
%% Simulate CO2 injection
sol = migrateInjection(Gt, res, pd, bestSingleCell, 'amount', 10, 'Nm', 1000, 'Ni', 50, 'T_migration', 2000*year);

%% Plot the actual migration
close all
target_trees = trees(arrayfun(@(x) any(traps == x.root), trees));

plotCellData(Gt, sol.h, sol.h > 0.1)
colors = {'r', 'y'};
for i = 1:2
    plotGrid(Gt, ismember(res.traps, target_trees(i).traps), 'FaceColor', colors{i}, 'EdgeColor', 'none');
end
plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', .05)
view(-210, 28)

%%
close all
plotGridBarchart(G, sol.h, 'edgecolor', 'none', 'facealpha', .6);
plotGrid(G, 'FaceAlpha', 0, 'edgealpha', .1)
fastRotateButton