%% Load results

[ws, states, reports] = getMultiplePackedSimulatorOutputs(problems);

%% Common parameters for plotting

% Plot wells
pw = @(G,W) plot(G.cells.centroids([W.cells], 1)     , ...
                 G.cells.centroids([W.cells], 2)     , ...
                 'ok', 'markerSize', 8, 'markerFaceColor', 'w', 'lineWidth', 2);

% Figures 
pos  = [-1000, 0, 500, 500];
posv = [-1000, 0, 500, 500];
fontSize = 12;
pth = fullfile(mrstPath('dg'), 'examples', 'rsc-2019', 'qfs-reorder', 'fig');
if 0
    savepng = @(name) print(fullfile(pth, name), '-dpng', '-r300');
    saveeps = @(name) print(fullfile(pth, name), '-depsc');
else
    savepng = @(name) [];
    saveeps = @(name) [];
end

hpos = [0.1300 0.1146 0.7750 0.0727];
cpos = [0.1300 0.07 0.7750 0.03];

% colors
gr = [1,1,1]*0.8;
clr = lines(2);

%% Plot saturation on refined grids

close all
frac = 0.8;
cmap = winter*frac + (1-frac);
timeSteps = [5, 15, 40, 60];
refFactor = [];
for tNo = timeSteps
    figure('Position', pos);
    hold on
    plotCellData(GF, states{1}{tNo}.order, 'edgec', 'none');
    pw(GF, WF)
    hold off
    axis equal tight; box on;
    ax = gca;
    [ax.XTickLabel, ax.YTickLabel] = deal({});
    colormap(jet)
    savepng(['qfs-reorder-', num2str(tNo)]);
end

%% Plot well solutions

close all
figure('Position', pos);
t = cumsum(schedule.step.val/day);
wcut = cellfun(@(ws) cellfun(@(ws) ws(2).wcut, ws), ws, 'unif', false);

hold on
plot(t, wcut{1}, 'linew', 2, 'color', clr(1,:));
plot(t, wcut{2}, '--', 'linew', 4, 'color', clr(2,:));
box on
axis([0, t(end), 0, 1])
xlabel('Time (days)');
ylabel('Water cut');
ax = gca;
ax.FontSize = fontSize;
legend({'Adaptive', 'Reference'}, 'location', 'northwest');
saveeps('qfs-adapt-wcut');

%%