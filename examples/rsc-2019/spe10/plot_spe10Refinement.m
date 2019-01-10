%% Load results

[ws, states, reports] = deal(cell(numel(degree), 4));

mdlIx = 1:4;
for dNo = 1:numel(degree)
    for mNo = mdlIx
        [ws{dNo, mNo}, states{dNo, mNo}, reports{dNo, mNo}] ...
            = getPackedSimulatorOutput(problems{dNo, mNo}, 'readFromDisk', false);
    end
end


%% Get iterations

iterations = cell(numel(degree), 4);
for dNo = 1:numel(degree)
    for mNo = 1:3
        if contains(names{mNo}, 'reorder') && ~isempty(reports{dNo, mNo})
            n   = reports{dNo,mNo}.numelData;
            rep = cell(n,1);
            for sNo = 1:n
                rep{sNo} = reports{dNo,mNo}{sNo};
            end
            iterations{dNo, mNo} = getReorderingTransportIterations(rep);
            reports{dNo, mNo} = rep;
        end
    end
end 

%% Common parameters for plotting

% Plot wells
pw = @(G,W) plot3(G.cells.centroids([W.cells], 1)     , ...
                  G.cells.centroids([W.cells], 2)     , ...
                  G.cells.centroids([W.cells], 3) + -3, ...
                 'ok', 'markerSize', 8, 'markerFaceColor', 'w', 'lineWidth', 2);

% pwn = @(G,W) text(G.cells.centroids([W(1).cells], 1) + 10     , ...
%                   G.cells.centroids([W(1).cells], 2) + 10    , ...
%                   G.cells.centroids([W(1).cells], 3) + -3, 'Inj', 'fontSize', 14);


% Figures 
pos  = [-1000, 0, 800, 500];
posv = [-1000, 0, 500, 800];
fontSize = 12;
pth = fullfile(mrstPath('dg'), 'examples', 'rsc-2019', 'spe10', 'fig');
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

%% Plot permeability

close all

figure('name', 'Perm', 'position', pos);

plotGrid(GC, 'facec', 'none');

perm = rock.perm(:,1);
plotCellData(G, log10(perm), 'edgec', 'none');
axis equal tight
view([90, 90])
colormap(pink)
hold on; pw(G, WF);

ax = gca;
ax.FontSize = fontSize;

[c, h] = colorbarHist(perm, [min(perm), max(perm)], 'South', 50, true);
h.Position = hpos;
c.Position = cpos;
c.Ticks = linspace(min(log10(perm)), max(log10(perm)), 5);
lbl = (10.^c.Ticks)/(100*milli*darcy);
c.TickLabels = round(lbl*1000)/1000;
c.FontSize = fontSize;

savepng('spe10-perm');

%% Plot porosity

close all

figure('name', 'Perm', 'position', pos);

plotGrid(GC, 'facec', 'none');

poro = rock.poro;
plotCellData(G, poro, 'edgec', 'none');
axis equal tight
view([90, 90])
colormap(pink)
hold on; pw(G, WF);

ax = gca;
ax.FontSize = fontSize;

[c, h] = colorbarHist(poro, [min(poro), max(poro)], 'South', 50);
h.Position = hpos;
c.Position = cpos;
c.Ticks = linspace(min(poro), max(poro), 5);
lbl = c.Ticks;
c.TickLabels = round(lbl*1000)/1000;
c.FontSize = fontSize;

savepng('spe10-poro');

%% Plot Saturation profiles for dG(0) and dG(1) reordered
cmap = mrstColormap();

close all

rIx = find(contains(names, 'reorder')); rIx = rIx(1);
st  = states(:, rIx);
rep = reports(:, rIx);
its = iterations(:, rIx);

timeSteps = [20, 42, 100];
itMaxdG0 = max(cellfun(@(it) max(it), its{1}(timeSteps)));
itMaxdG1 = max(cellfun(@(it) max(it), its{2}(timeSteps)));
itMax = max(itMaxdG0, itMaxdG1);
% itMax = cellfun(@(it) max(it), cellfun(@(it) max(cellfun(@(it) max(it(timeSteps)), it)), its, 'unif', false));
% itMax = max(itMax(1:2));
cmap = winter;
cmap = cmap(end:-1:1, :);
for tNo = timeSteps
    for dNo = 1:2%numel(degree)
        if ~isempty(st{dNo})

            figure('position', posv, 'name', ['dG(', num2str(degree(dNo)), ')']);

%             plotGrid(G, 'facec', 'none', 'edgec', 'none')
%             c = its{dNo}{tNo} > 0 | st{dNo}{tNo}.s(:,1) > 0.21;
%             plotCellData(G, st{dNo}{tNo}.s(c,1), c, 'edgec', 'none'); 
            unstructuredContour(G, st{dNo}{tNo}.s(:,1), 10,'linew', 2);
            hold on
            pw(G, WF);
            axis equal tight
            box on
            caxis([0.2, 0.8]);
            ax = gca;
            [ax.XTickLabel, ax.YTickLabel] = deal({});
            colormap(cmap)
            savepng(['spe10-sat-', num2str(tNo), '-dg', num2str(degree(dNo))]);
            
            figure('position', posv, 'name', ['dG(', num2str(degree(dNo)), ') iterations']);

            plotGrid(G, 'facec', 'none', 'edgec', 'none')
            c = its{dNo}{tNo} > 0;
            plotCellData(G, its{dNo}{tNo}(c), c, 'edgec', 'none'); 
            hold on
            pw(G, WF);
            axis equal tight
            box on
            caxis([0, itMax]);
            ax = gca;
            [ax.XTickLabel, ax.YTickLabel] = deal({});
            colormap(jet)
            savepng(['spe10-its-', num2str(tNo), '-dg', num2str(degree(dNo))]);
            

        end
    end
end

%% Plot adaptive refinement for dG(0) and dG(1)

close all

aIx = find(contains(names, 'adapt')); aIx = aIx(1);
st  = states(:, aIx);
rep = reports(:, aIx);
its = iterations(:, aIx);
timeSteps = [20, 42, 60];

for tNo = timeSteps
    for dNo = 1:2
        if ~isempty(st{dNo})

            figure('position', posv, 'name', ['dG(', num2str(degree(dNo)), ')']);

%              pink = [214, 154, 153]/255;
            hold on
            plotGrid(st{dNo}{tNo}.G, 'facec', 'none', 'edgecolor', gray);
            unstructuredContour(G, st{dNo}{tNo}.s(:,1), 10, 'linewidth', 2);
%             plotCellData(G, st{dNo}{tNo}.s(:,1), 'edgec', 'none'); 
            pw(G, WF);
            axis equal tight
            box on
            caxis([0.2, 0.8]);
            ax = gca;
            [ax.XTickLabel, ax.YTickLabel] = deal({});
            colormap(cmap)
            savepng(['spe10-ref-', num2str(tNo), '-dg', num2str(degree(dNo))]);

            

        end
    end
end
