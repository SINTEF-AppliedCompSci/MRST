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
timeSteps = [10, 40, 60, 80, 108];
stRef    = states{2};
stAdapt  = states{4};
stCoarse = states{5};

% its = iterations{2};
% itMax = max(cellfun(@max, its));
G  = computeCellDimensions2(model.G);
GC = generateCoarseGrid(G, coarsemodel.transportModel.G.partition);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC);


for tNo = timeSteps
    
%     c  = it > 0;

    if ~isempty(stRef{tNo})
        figure('position', pos, 'Name', 'Reference')
        hold on
        unstructuredContour(model.G, stRef{tNo}.x(:,1), 5, 'linew', 2);
        colormap(copper);
        axis equal tight
        pw(G, W)
        ax = gca;
        [ax.XTickLabel, ax.YTickLabel] = deal({});
        box on
        drawnow();pause(0.5)
        savepng(['qfs-co2-x-', num2str(tNo), '-ref']);
    end
    
    if ~isempty(stCoarse{tNo})
        figure('position', pos, 'Name', 'Coarse')
        hold on
        unstructuredContour(GC, stCoarse{tNo}.x(:,1), 5, 'linew', 2);
        colormap(copper);
        axis equal tight
        pw(G, W)
        ax = gca;
        [ax.XTickLabel, ax.YTickLabel] = deal({});
        box on
        drawnow();pause(0.5)
        savepng(['qfs-co2-x-', num2str(tNo), '-coarse']);
    end
    
    if ~isempty(stAdapt{tNo})
        figure('position', pos, 'Name', 'Adaptive')
        hold on
        unstructuredContour(model.G, stAdapt{tNo}.x(:,1), 5, 'linew', 2);
        plotGrid(stAdapt{tNo}.G, 'facec', 'none', 'edgec', gray);
        colormap(copper);
        axis equal tight
        pw(G, W)
        ax = gca;
        [ax.XTickLabel, ax.YTickLabel] = deal({});
        box on
        drawnow();pause(0.5)
        savepng(['qfs-co2-x-', num2str(tNo), '-adapt']);
    end
            
end

disp(realTime(timeSteps));

%% Well curves

close all

comp = cell(6,1);
for mNo = 1:6
    cmp = cellfun(@(ws) ws(2).components, ws{mNo}, 'unif', false);
    comp{mNo} = abs(vertcat(cmp{:}));
end

lt = {'', '-', '--', '.', '.'};
lw = [2,2,4,2,2];
mz = 10*ones(5,1);
mz(end) = 10;
names = cellfun(@(p) p.Description, problems, 'unif', false);

dtt = cumsum(schedule.step.val)/day;
clr = copper(3);
for cNo = 1:size(comp{2},2)
    figure();
    hold on
    cMax = -Inf;
    for mNo = 2:6
        cmp = comp{mNo};
        if ~isempty(cmp)
            cMax = max(cMax, max(cmp(:,cNo)));
            plot(dtt(3:size(cmp,1)), cmp(3:end,cNo), lt{mNo}, 'linew', lw(mNo), 'markerSize', mz(mNo));
        end
    end
    hold off
    legend(names(2:6), 'location', 'northwest');
    box on; grid on;
    ax = gca;
    ax.FontSize = fontSize;
    axis([0, dtt(end), 0, cMax]);
    xlabel('Time (days)')
    ylabel('Rate (m^3/s)')
    saveeps(['qfs-co2-well-', num2str(cNo)]);
    
end

%%

iterations = cell(numel(problems),1);

for pNo = 2:numel(problems)
    
    if ~isempty(reports{pNo})
        if contains(problems{pNo}.Description, 'Reordering')
            iterations{pNo} = getReorderingTransportIterations(reports{pNo});
        else
            for sNo = 1:numel(reports{pNo})
                if isfield(states{pNo}{sNo}, 'G')
                    g = states{pNo}{sNo}.G;
                elseif pNo == 5
                    g = GC;
                else
                    g = model.G;
                end
                iterations{pNo}{sNo} = ...
                    repmat(reports{pNo}{sNo}.StepReports{1}.NonlinearReport{1}.TransportSolver.Iterations, g.cells.num, 1);
            end
        end
    end
end

intIts = cellfun(@(its) cellfun(@sum, its), iterations(2:end-1), 'unif', false);

figure('position', [-1000, 0, 600, 300]);
hold on
for pNo = 1:4
    plot(dtt, intIts{pNo}, '-', 'linew', 2)
end
hold off
box on; grid on
xlim([0, dtt(end)])
legend({'Sequential', 'Reordered', 'Adaptive', 'Coarse'});
ax = gca;
ax.FontSize = fontSize;
xlabel('Time (days)');
ylabel('Transport iterations');
saveeps('qfs-co2-iterations');

