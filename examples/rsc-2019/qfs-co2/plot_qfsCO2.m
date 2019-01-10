%% Load results

[ws, states, reports] = getMultiplePackedSimulatorOutputs(problems, 'readFromDisk', false);

%% Get iterations

rIx = 3;
iterations = cell(3,1);
for mNo = 3
    n   = reports{mNo}.numelData;
    rep = cell(n,1);
    for sNo = 1:n
        rep{sNo} = reports{mNo}{sNo};
    end
    iterations{mNo} = getReorderingTransportIterations(rep);
    reports{mNo} = rep;
end

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

realTime = cumsum(schedule.step.val)/day;

close all

name = cellfun(@(p) p.Name, problems, 'unif', false);
st  = states{rIx};
timeSteps = [10, 40, 60, 80, 108];

its = iterations{3};
itMax = max(cellfun(@max, its));


for tNo = timeSteps
    
    it = its{tNo};
    c  = it > 0;
    
    figure('position', pos)
%     unstructuredContour(model.G, st{tNo}.x(:,1), 20, 'linew', 2);
    plotCellData(model.G, st{tNo}.x(:,1), 'edgec', 'none');
    colormap(copper);
    axis equal tight
    hold on
    pw(G, W)
    ax = gca;
    [ax.XTickLabel, ax.YTickLabel] = deal({});
    box on
    savepng(['qfs-co2-x-', num2str(tNo)]);
    
    figure('position', pos)
    plotCellData(model.G, it(c), c, 'edgec', 'none');
    colormap(jet);
    axis equal tight
    hold on
    pw(G, W)
    caxis([1,itMax]);
    ax = gca;
    [ax.XTickLabel, ax.YTickLabel] = deal({});
    box on
    savepng(['qfs-co2-its-', num2str(tNo)]);
    
end

disp(realTime(timeSteps));

%%

close all
plotWellSols(ws(1:2))

%%

