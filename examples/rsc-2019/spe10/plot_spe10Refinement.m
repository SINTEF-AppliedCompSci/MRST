%% Load results

[ws, states, r] = deal(cell(numel(degree), 4));

mdlIx = 1:4;
for dNo = 1:numel(degree)
    for mNo = mdlIx
        [ws{dNo, mNo}, states{dNo, mNo}, r{dNo, mNo}] ...
            = getPackedSimulatorOutput(problems{dNo, mNo}, 'readFromDisk', false);
    end
end
[wsWENO, statesWENO, reportsWENO] = getPackedSimulatorOutput(weno, 'readFromDisk', false);

%% Get iterations

reports = cell(size(r));
iterations = cell(numel(degree), 4);
for dNo = 1:numel(degree)
    for mNo = 1:3
        n   = r{dNo,mNo}.numelData;
        rep = cell(n,1);
        for sNo = 1:n
            rep{sNo} = r{dNo,mNo}{sNo};
        end
        if contains(names{mNo}, 'reorder') && ~isempty(r{dNo, mNo})
            iterations{dNo, mNo} = getReorderingTransportIterations(rep);
        else
            for sNo = 1:numel(rep)
                st = states{dNo, mNo}{sNo};
                if isfield(st, G); g = st.G; else; g = G; end
                its = repmat(rep{sNo}.StepReports{1}.NonlinearReport{1}.TransportSolver.Iterations, g.cells.num, 1);
                iterations{dNo, mNo}{sNo} = its;
            end
        end
        reports{dNo, mNo} = rep;
    end
end 

%% Common parameters for plotting

% Plot wells
pw = @(G,W) plot3(G.cells.centroids([W.cells], 1)     , ...
                  G.cells.centroids([W.cells], 2)     , ...
                  G.cells.centroids([W.cells], 3) + -3, ...
                 'ok', 'markerSize', 8, 'markerFaceColor', 'w', 'lineWidth', 2);

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
clr = lines(5);

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

close all

rIx = find(contains(names, 'reorder')); rIx = rIx(1);
st  = states(:, 1);
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
            
            % Plot dG profile
            figure('position', posv, 'name', ['dG(', num2str(degree(dNo)), ')']);

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
            
            % Plot reordering iterations
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
    
    if ~isempty(statesWENO{tNo})
        % Plot WENO profile
        figure('position', posv, 'name', 'WENO');

        unstructuredContour(G, statesWENO{tNo}.s(:,1), 10, 'linew', 2);
        hold on
        pw(G, WF);
        axis equal tight
        box on
        caxis([0.2, 0.8]);
        ax = gca;
        [ax.XTickLabel, ax.YTickLabel] = deal({});
        colormap(cmap)
        savepng(['spe10-sat-', num2str(tNo), '-weno']);
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

%% Plot Iterations

close all
intIts = cellfun(@(its) cellfun(@sum, its), iterations(1:2,1:3), 'unif', false);
intItsWENO = cellfun(@(r) r.Iterations*G.cells.num, {reportsWENO{1:reportsWENO.numelData}});
dtt    = cumsum(schedule.step.val)/day;

itPos = [-1000, 0, 600, 300];
for dNo = 1:2
    figure('position', itPos)
    pNames = {['dG(', num2str(dNo-1), ')']           , ...
              ['dG(', num2str(dNo-1), '), reordered'], ...
              ['dG(', num2str(dNo-1), '), adaptive'] , ...
              'WENO'};
    hold on
    imax = -Inf;
    for mNo = 1:3
        imax = max(max(intIts{dNo,mNo}), imax);
        plot(dtt, intIts{dNo,mNo}, '-', 'linew', 2)
    end
    if dNo == 2
        plot(dtt, intItsWENO, '-', 'linew', 2)
    end
    hold off
    legend(pNames)
    box on; grid on
    axis([0, dtt(end), 0, imax*1.1]);
    ax = gca;
    ax.FontSize = fontSize;
    xlabel('Time (days)');
    ylabel('Iterations');
    saveeps(['spe10-iterations-', num2str(dNo-1)]);
end

%% Plot well curves

[wellSols, wcut] = deal(cell(2, 3));
for dNo = 1:2
    for mNo = 1:3
        if ws{dNo, mNo}.numelData > 0
            wellSols{dNo, mNo} =  {ws{dNo, mNo}{1:ws{dNo,mNo}.numelData}};
            wcut{dNo, mNo} = cellfun(@(ws) ws(2).wcut, wellSols{dNo,mNo});
        end
    end
end
wellSolsWENO = {wsWENO{1:wsWENO.numelData}};
wcutWENO     = cellfun(@(ws) ws(2).wcut, wellSolsWENO);

%%

close all
wcutPos = [-1000, 0, 400, 300];

mStyle = {'-', '--', '.'};
mSize  = 10;
lw = [2,4,2];
it = 1;
for dNo = 1:2
    figure('position', wcutPos)
    pNames = {['dG(', num2str(dNo-1), ')']           , ...
              ['dG(', num2str(dNo-1), '), reordered'], ...
              ['dG(', num2str(dNo-1), '), adaptive'] , ...
              'WENO'};
    hold on
    for mNo = 1:3
        plot(dtt, wcut{dNo,mNo}, mStyle{mNo}, 'linew', lw(mNo), 'markerSize', mSize);
    end
    hold off
    legend(pNames, 'location', 'northwest')
    box on; grid on
    axis([0, dtt(end), 0, 1]);
    ax = gca;
    ax.FontSize = fontSize;
    xlabel('Time (days)');
    ylabel('Water cut');
    saveeps(['spe10-wcut-dg', num2str(dNo-1)]);
end


figure('position', wcutPos)
hold on
for dNo = 1:2
    for mNo = 1
        plot(dtt, wcut{dNo,mNo}, mStyle{mNo}, 'linew', lw(mNo), 'markerSize', mSize);
    end
end
plot(dtt, wcutWENO, 'linew', lw(mNo), 'markerSize', mSize, 'color', clr(4,:));
hold off
box on; grid on
axis([0, dtt(end), 0, 1]);
legend({'dG(0)', 'dG(1)', 'WENO'}, 'location', 'northwest');
ax = gca;
ax.FontSize = fontSize;
xlabel('Time (days)');
ylabel('Water cut');
saveeps('spe10-wcut');

%%


























