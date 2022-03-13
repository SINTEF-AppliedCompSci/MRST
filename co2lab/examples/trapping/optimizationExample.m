%% Optimal injection points
% In this example we demonstrate how information about structural traps and
% their connections can be used to define a good guess of injector
% placements. That is, we will consider the problem of determining where N
% wells should be placed to optimize the potential for structural trapping.
mrstModule add co2lab;

%% Extract a subset of the Utsira formation for analysis
coarsefactor = 2;
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

res = trapAnalysis(Gt, true);
trees = maximizeTrapping(Gt, 'res', res, 'removeOverlap', true);
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
mrstModule add mrst-gui;
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
    plotCellData(G, repmat(v(i), G.cells.num, 1), cells & ~plotted, ...
                 'EdgeColor','none');
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
camlight head
lighting phong
% set(gca,'ZDir','reverse');

axis tight off
view(64, 72)
caxis([min(treestart), max(treestart)]);

%% Alternative visualization
% This is more dull, but has a somewhat higher information content
n = max(res.trap_regions);
c = jet(n); c=c(randperm(n),:);

subplot(1,2,2);
plotted = false(G.cells.num, 1);
for i = 1:numel(trees)
    cells = ismember(res.traps, trees(i).traps);
    plotCellData(G, repmat(v(i), G.cells.num, 1), cells & ~plotted, ...
                 'EdgeColor','none');
    plotted = plotted | cells;
end
colormap(jet(128));
contourAtlas(d{2},50,.25,.4*ones(1,3));
coord = Gt.cells.centroids(trapcells,:);
plot(coord(:,1),coord(:,2),'o',...
   'MarkerEdgeColor','b','MarkerFaceColor','r','LineWidth',3,'MarkerSize',8);

perc = @(x) 100*(10^x)/sum(tv);
plotBar(treestart, trapcells, @(x) sprintf('%1.1f%%', perc(x)))


subplot(1,2,1);
for i=1:n
   plotFaces(G,boundaryFaces(G,res.trap_regions==i), ...
      'EdgeColor','none','FaceColor', .5*c(i,:)+.5*ones(1,3));
end
plotFaces(G,boundaryFaces(G,res.trap_regions==0), ...
   'EdgeColor','none','FaceColor', .85*ones(1,3));
plotted = false(G.cells.num, 1);
for i = 1:numel(trees)
    cells = ismember(res.traps, trees(i).traps);
    plotCellData(G, repmat(v(i), G.cells.num, 1), cells & ~plotted, ...
                 'EdgeColor', 'none');
    plotted = plotted | cells;
end
colormap(jet(128));

%%
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
