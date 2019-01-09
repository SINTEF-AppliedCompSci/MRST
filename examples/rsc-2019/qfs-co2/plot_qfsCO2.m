%% Load results

[ws, states, reports] = getMultiplePackedSimulatorOutputs(problems);

%% Get iterations

rIx = strcmpi(name, 'reorder');
iterations = getReorderingTransportIterations(reports{3});

%%

G = model.G;
W = problems{1}.SimulatorSetup.schedule.control(1).W;
t = cumsum(problems{1}.SimulatorSetup.schedule.step.val);

% Plot wells
pw = @(G,W) plot3(G.cells.centroids([W.cells], 1)     , ...
                  G.cells.centroids([W.cells], 2)     , ...
                  G.cells.centroids([W.cells], 3) + -10, ...
                 'ok', 'markerSize', 8, 'markerFaceColor', 'w', 'lineWidth', 2);

% Figures 
pos = [-1000, 0, 500, 500];
fontSize = 12;
pth = fullfile(mrstPath('dg'), 'examples', 'rsc-2019', 'qfs-co2', 'fig');
if 1
    savepng = @(name) print(fullfile(pth, name), '-dpng', '-r300');
    saveeps = @(name) print(fullfile(pth, name), '-depsc');
else
    savepng = @(name) [];
    saveeps = @(name) [];
end

hpos = [0.1300 0.1146 0.7750 0.0727];
cpos = [0.1300 0.07 0.7750 0.03];

%% Plot compositions and 

close all

name = cellfun(@(p) p.Name, problems, 'unif', false);
st  = states{rIx};
timeSteps = [5, 10, 25];

for tNo = timeSteps
    
    it = iterations{tNo};
    c  = it > 0;
    
    figure('position', pos)
    plotCellData(model.G, st{tNo}.x(c,2), c);
    colormap(bone);
    caxis([0,1]);
    axis equal tight
    hold on
    pw(G, W)
    ax = gca;
    [ax.XTickLabel, ax.YTickLabel] = deal({});
    box on
    savepng(['qfs-co2-x-', num2str(tNo)]);
    
    figure('position', pos)
    plotCellData(model.G, it(c), c);
    colormap(jet);
    axis equal tight
    hold on
    pw(G, W)
    ax = gca;
    [ax.XTickLabel, ax.YTickLabel] = deal({});
    box on
    savepng(['qfs-co2-its-', num2str(tNo)]);
    
end

%%

close all
plotWellSols(ws(1:2))

%%

