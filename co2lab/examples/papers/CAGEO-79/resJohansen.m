%% Resolution effects: Johansen formation
% When working with subsea reservoirs coarsening will always be a factor:
%
% * Simulations on very large grids become prohibitively expensive in terms
% of computing power very fast.
% * All grids are essentially coarse approximations of highly complex
% geometry. Different grids have a different scale for each cell depending on
% what kind of data was used to produce the grid in the first place and the
% intended usage of the datasets.
%
% This example demonstrates the effects of geometry coarsening on the CO2
% Storage Atlas grids, with a special focus on structural trapping. To gain a
% rough of estimate the information missing in the coarse Atlas grids, we will
% also compare the relatively fine scale Sleipner grid with the same region in
% the Utsira formation, where Sleipner is located.
mrstModule('add', 'coarsegrid', 'deckformat', 'mex', 'matlab_bgl', 'opm_gridprocessing', ...
            'libgeometry');
   

%% Generate grids and trapping for the Johansen formation
% We generate six realizations of the Johansen formation from the CO2 Storage
% Atlas: The first is the full dataset, the second coarsened by a factor 2 in
% both i and j directions, the third a factor 3, and so on.  This gives a set of
% grids where the most detailed grid has approximately 80,000 cells while the
% smallest has ~2,000 cells. In the finest grid each cell has a resolution of
% 500x500 m^2 per cell whereas the coarsest has 3000x3000 m^2: Both are fairly
% large in terms of typical simulation grids, but the size of the cells in the
% coarsest grid is several orders of magnitude larger than in the finest.

N = 6;
[Grids, res] = deal(cell(numel(N),1));

tic;
for i = 1:N
    gr = getAtlasGrid('Johansenfm', 'coarsening', i);
    % Create top surface grid
    try
       G = processgrid(gr{1});
       Gt = topSurfaceGrid(mcomputeGeometry(G(1)));
    catch
       G = processGRDECL(gr{1});
       Gt = topSurfaceGrid(computeGeometry(G(1)));
    end
    Grids{i} = Gt;
    % Create trapping and store volumes of each trap
    res{i} = trapAnalysis(Gt, true);
    res{i}.volumes = volumesOfTraps(Gt, res{i}, []);
end
toc;

%% Number and size distribution of global traps
% Output number of global traps, as well as their average volumes, for each
% degree of coarsening.  Also plot their size distribution.
for i = 1:N
    fprintf('Coarsening level %d:\n', i);
    fprintf('  Num. global traps: %d\n', numel(res{i}.volumes));
    fprintf('  Total trap volume: %e m3\n', sum(res{i}.volumes));
    fprintf('  Avg. global trap size: %e m3\n', mean(res{i}.volumes));
end

figure;
defaultpos = get(0, 'DefaultFigurePosition');
set(gcf, 'Position', [defaultpos(1:2) - [300 100], 1200 800],...
   'PaperPositionMode','auto');
colorize = 'rgcmyb';

% Plot the total volume as subplot
subplot(2,2,1); cla
hold on
vol = cellfun(@(x) sum(x.volumes), res);
for i = 1:N
    bar(i, vol(i), colorize(i))
end
title('Total trap volume [m^3]', 'FontSize', 14)
set(gca, 'Color', get(gcf, 'Color'), 'FontSize', 14, ...
   'XTickLabel', regexp(num2str((1:N)*500), '  ', 'split'))
xlabel('Lateral resolution [m]', 'FontSize', 12);
axis tight

% Plot cumulative colume covered
subplot(2, 2, 2); cla
hold on;
cum_max = 0;
for i = 1:N
    cumul_vol = cumsum(sort(res{i}.volumes, 'descend'));
    plot([0, cumul_vol], colorize(i),'LineWidth', 1, 'Marker', '.', 'MarkerSize', 16);
    cum_max = max(cum_max, cumul_vol(end));
end
axis([0 50 0 cum_max*1.05]);
set(gca, 'Color', get(gcf,'Color'), 'FontSize', 14);
h=legend(regexp(num2str((1:N)*500), '  ', 'split'), 'location', 'East');
title('Cumulative trap volume [m^3]', 'FontSize', 14);
xlabel('Number of traps counted', 'FontSize', 12);

% Bar plot: number of traps at each coarsening level
subplot(2, 2, 3);
hold on
for i = 1:N
    bar(i, numel(res{i}.volumes), colorize(i));
end
title('Total number of traps', 'FontSize', 14)
axis tight;
set(gca, 'Color', get(gcf, 'Color'), 'FontSize', 14, ...
   'XTick', [1 2 3 4 5 6], ...
   'XTickLabel', regexp(num2str((1:N)*500), '  ', 'split'))
xlabel('Lateral resolution [m]', 'FontSize', 12);

% Plot the average volume as subplot
subplot(2,2,4); cla
hold on
mvol = cellfun(@(x) mean(x.volumes), res);
for i = 1:N
    bar(i, mvol(i), colorize(i))
end
title('Average trap volume [m^3]', 'FontSize', 14)
set(gca, 'Color', get(gcf, 'Color'), 'FontSize', 14, ...
   'XTickLabel', regexp(num2str((1:N)*500), '  ', 'split'))
axis tight
xlabel('Lateral resolution [m]', 'FontSize', 12);

%% Plot a subset of Johansen with different degree of coarsening
% We first define a subdomain consisting of a minimum and maximum value for
% both x and y coordinates which is then plotted on the fine grid.
G = Grids{1};

subdomain = [5.25e5, 6.70e6; 5.30e5, 6.75e6];

x = G.cells.centroids(:,1);
y = G.cells.centroids(:,2);

subset = x > subdomain(1,1) & x < subdomain(2,1) &...
         y > subdomain(1,2) & y < subdomain(2,2);
figure; clf;
plotCellData(G, G.cells.z, 'EdgeColor','none')
[~,dset] = getAtlasGrid('Johansenfm', 'coarsening', 1);
contourAtlas(dset{2},20,.5,'k');
plotGrid(G, subset, 'EdgeColor', 'None', 'FaceColor', 'black', 'FaceAlpha', .6)
axis tight, box on, set(gca,'Color','none'); 
colorbar('Location','EastOutside','FontSize',16);

%% The effect of coarsening on trapping analysis
% It is obvious that fine structural details are lost when coarsening the
% grids. The coarsening operation acts as a smoother on the grid, removing
% ridges, folds and oscillations that are present on a shorter wavelength
% than the coarse cells. Unfortunately, these small oscillations are
% especially interesting for CO2 migration studies: Small local height
% maxima can divert small "rivers" of CO2 and act as structural traps.
%
% To demonstrate that these traps are lost when coarsening, we plot the
% structural traps estimated by trapAnalysis for the different grids. Note
% that several smaller traps are removed as the coarsening increases, which
% can be shown statistically by noting that the mean of the trap volume
% quickly increases as the smaller traps are smoothed away. 
%
% The total trapping volume also changes as the coarsening is increased: In
% the beginning the volume increases as the largest traps become slightly
% larger due to the lower resolution. Later on, the total volume shrinks as
% smaller traps are removed entirely.

% Plot the traps
clf
defaultpos = get(0, 'DefaultFigurePosition');
set(gcf, 'Position', [defaultpos(1:2) - [0 100], 800 400],...
   'PaperPositionMode','auto');
hl = [];
for i = 1:N
    G = Grids{i};
    G.nodes.z = G.nodes.z + 1000*i;
    
    x = G.cells.centroids(:,1);
    y = G.cells.centroids(:,2);
    subset = x > subdomain(1,1) & x < subdomain(2,1) &...
             y > subdomain(1,2) & y < subdomain(2,2);
         
    h = plotGrid(G, subset, 'facec', colorize(i), 'facea', .2, ...
                 'edgec', colorize(i), 'edgea', .9);
    hl = [hl; h];%#ok
    tr = res{i};
    G_flat = flattenTraps(G, tr);
    G_flat.nodes.z = G_flat.nodes.z + 1000*i;
    plotGrid(G_flat, subset & tr.traps ~= 0, 'FaceColor', colorize(i), 'EdgeColor', 'none')
end
% view(-30, 60)
view(86, 12)
axis tight off
title('Traps for successively coarsed grids','FontSize',14)
legend(hl, regexp(num2str((1:N)*500), '  ', 'split'), ...
   'Location','SouthOutside', 'Orientation', 'horiz');
legend boxoff


