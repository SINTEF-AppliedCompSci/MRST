%% Load results

[ws, states, reports] = getMultiplePackedSimulatorOutputs(problems);

%%

% Plot wells
pw = @(G,W) plot(G.cells.centroids([W.cells], 1)     , ...
                 G.cells.centroids([W.cells], 2)     , ...
                 'ok', 'markerSize', 8, 'markerFaceColor', 'w', 'lineWidth', 2);

% Figures 
pos  = [-1000, 0, 500, 500];
posv = [-1000, 0, 500, 500];
fontSize = 12;
pth = fullfile(mrstPath('dg'), 'examples', 'rsc-2019', 'qfs-adapt-reorder', 'fig');
if 1
    savepng = @(name) print(fullfile(pth, name), '-dpng', '-r300');
    saveeps = @(name) print(fullfile(pth, name), '-depsc');
else
    savepng = @(name) [];
    saveeps = @(name) [];
end

hpos = [0.1300 0.1146 0.7750 0.0727];
cpos = [0.1300 0.07 0.7750 0.03];

% colors
gray = [1,1,1]*0.5;
clr = lines(2);

cmap = winter;
cmap = cmap(end:-1:1,:);

%%

figure('Position', pos);
hold on
unstructuredContour(state.G, state.s(:,1), 5, 'linew', 2);
plotGrid(state.G, 'facec', 'none', 'edgec', gray);
caxis([0,1]);
pw(state.G, forces.W)
hold off
axis equal tight; box on;
ax = gca;
[ax.XTickLabel, ax.YTickLabel] = deal({});
colormap(cmap)
% savepng(['qfs-adapt-', num2str(tNo)]);
% refFactor = [refFactor, rf(states{1}{tNo}.G)];