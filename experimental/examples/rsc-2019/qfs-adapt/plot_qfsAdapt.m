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
fontSize = 14;
pth = fullfile(mrstPath('dg'), 'examples', 'rsc-2019', 'qfs-adapt', 'fig');
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
clr = lines(3);

%% Plot saturation on refined grids

close all
rf = @(G) G.cells.num./GF.cells.num;
frac = 0.8;
cmap = winter*frac + (1-frac);
timeSteps = [5, 15, 30, 50];
refFactor = [];
for tNo = timeSteps
    figure('Position', pos);
    hold on
    unstructuredContour(GF, states{1}{tNo}.s(:,1), 5, 'linew', 2);
%     plotCellData(GF, states{1}{tNo}.s(:,2), 'edgec', 'none');
    plotGrid(states{1}{tNo}.G, 'facec', 'none', 'edgec', gray);
    caxis([0,1]);
    pw(GF, WF)
    hold off
    axis equal tight; box on;
    ax = gca;
    [ax.XTickLabel, ax.YTickLabel] = deal({});
    colormap(cmap)
    savepng(['qfs-adapt-', num2str(tNo)]);
    refFactor = [refFactor, rf(states{1}{tNo}.G)];
end

%% Plot well solutions

close all
figure('Position', pos);
t = cumsum(schedule.step.val/day);
wcut = cellfun(@(ws) cellfun(@(ws) ws(2).wcut, ws), ws, 'unif', false);

hold on
plot(t, wcut{2}, 'linew', 2, 'color', clr(1,:));
plot(t, wcut{1}, '--', 'linew', 4, 'color', clr(2,:));
plot(t, wcut{3}, 'linew', 3, 'color', clr(3,:));
box on
axis([0, t(end), 0, 1])
xlabel('Time (days)');
ylabel('Water cut');
ax = gca;
ax.FontSize = fontSize;
legend({'Reference', 'Adaptive', 'Coarse'}, 'location', 'northwest');
saveeps('qfs-adapt-wcut');

%%