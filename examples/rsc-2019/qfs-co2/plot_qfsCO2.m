%% Load results

[ws, states, reports] = getMultiplePackedSimulatorOutputs(problems, 'readFromDisk', true);

%% Get iterations

rIx = 3;
iterations = cell(4,1);
for mNo = 1:4
%     n   = reports{mNo}.numelData;
    rep = reports{mNo};
%     rep = cell(n,1);
%     for sNo = 1:n
%         rep{sNo} = reports{mNo}{sNo};
%     end
    if mNo == 3
        iterations{mNo} = getReorderingTransportIterations(rep);
    elseif mNo > 1
        iterations{mNo} = cellfun(@(rep) rep.StepReports{1}.NonlinearReport{1}.TransportSolver.Iterations, rep);
    end
    
%     reports{mNo} = rep;
end

%% Common

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

gray = [1,1,1]*0.5;

%% Plot co2 compositions for adaptive

realTime = cumsum(schedule.step.val)/day;

close all

name = cellfun(@(p) p.Name, problems, 'unif', false);
st  = states{4};
timeSteps = [10, 40, 60, 80, 108];
stref = states{2};
stc   = states{5};

% its = iterations{2};
% itMax = max(cellfun(@max, its));
G  = computeCellDimensions2(model.G);
GC = generateCoarseGrid(G, coarsemodel.transportModel.G.partition);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC);


for tNo = timeSteps
    
    c  = it > 0;
    
    figure('position', pos)
    hold on
    unstructuredContour(model.G, st{tNo}.x(:,1), 5, 'linew', 2);
    plotGrid(st{tNo}.G, 'facec', 'none', 'edgec', gray);
    colormap(copper);
    axis equal tight
    pw(G, W)
    ax = gca;
    [ax.XTickLabel, ax.YTickLabel] = deal({});
    box on
    drawnow();pause(0.5)
    savepng(['qfs-co2-x-', num2str(tNo), '-adapt']);
    
    figure('position', pos)
    hold on
    unstructuredContour(model.G, stref{tNo}.x(:,1), 5, 'linew', 2);
    colormap(copper);
    axis equal tight
    pw(G, W)
    ax = gca;
    [ax.XTickLabel, ax.YTickLabel] = deal({});
    box on
    drawnow();pause(0.5)
    savepng(['qfs-co2-x-', num2str(tNo), '-ref']);
    
    figure('position', pos)
    hold on
    unstructuredContour(GC, stc{tNo}.x(:,1), 5, 'linew', 2);
    colormap(copper);
    axis equal tight
    pw(G, W)
    ax = gca;
    [ax.XTickLabel, ax.YTickLabel] = deal({});
    box on
    drawnow();pause(0.5)
    savepng(['qfs-co2-x-', num2str(tNo), '-coarse']);
    
    
end

disp(realTime(timeSteps));

%% Iteration plot

%% Well curves

close all

comp = cell(6,1);
for mNo = 1:6
    cmp = cellfun(@(ws) ws(2).components, ws{mNo}, 'unif', false);
    comp{mNo} = abs(vertcat(cmp{:}));
end

lt = {'-', '-', '-.', '--', ':', '-o'};
names = cellfun(@(p) p.Description, problems, 'unif', false);

dtt = cumsum(schedule.step.val);
clr = copper(3);
for cNo = 1:size(comp{1},2)
    figure();
    hold on
    for mNo = 2:6
        cmp = comp{mNo};
        plot(dtt(3:size(cmp,1)), cmp(3:end,cNo), lt{mNo});
    end
    hold off
    legend(names(2:6));
end

%%

close all
plotWellSols(ws(2:6))

%%

% figure('position', pos)
%     plotCellData(model.G, it(c), c, 'edgec', 'none');
%     colormap(jet);
%     axis equal tight
%     hold on
%     pw(G, W)
%     caxis([1,itMax]);
%     ax = gca;
%     [ax.XTickLabel, ax.YTickLabel] = deal({});
%     box on
%     savepng(['qfs-co2-its-', num2str(tNo)]);

