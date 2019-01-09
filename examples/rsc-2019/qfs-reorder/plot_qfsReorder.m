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
timeSteps = [5, 15, 30, 40];
refFactor = [];
for tNo = timeSteps
    figure('Position', pos);
    hold on
    st = states{1}{tNo}; 
    plotCellData(G, states{1}{tNo}.order, 'edgec', 'none');
    pw(G, W)
    s = states{1}{tNo}.s(:,1);
    unstructuredContour(G, s, 'color', 'w', 'linew', 2);
    hold off
    axis equal tight; box on;
    ax = gca;
    [ax.XTickLabel, ax.YTickLabel] = deal({});
    colormap(jet)
    savepng(['qfs-reorder-', num2str(tNo)]);
    
%     figure('Position', pos);
    order = st.order;
    c = sum(order == order',2) > 1;
%     plotGrid(G, 'facec', 'none');
    plotGrid(G, c, 'facec', 'white');
   
end

%%