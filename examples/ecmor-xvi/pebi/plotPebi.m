%% Load data

[statesDG, statesDGReorder, wcutDG, wcutDGReorder, tDG, tDGReorder, itDG, itDGReorderAll, itDGReorder] = deal(cell(numel(degree),1));

useH = false;
ext = [];
if useH
    ext = '-homogeneous';
end

for dNo = 1:2%numel(degree)
    
    ohDG = getOutHandler(['dg', num2str(degree(dNo)), ext]);
    ohDGReorder = getOutHandler(['dg', num2str(degree(dNo)), '-reorder', ext]);
    
    rhDG = getRepHandler(['dg', num2str(degree(dNo)), ext]);
    rhDGReorder = getRepHandler(['dg', num2str(degree(dNo)), '-reorder', ext]);

    for sNo = 1:25%ohDG.numelData

        st = ohDG{sNo};
        statesDG{dNo}{sNo} = st;
        wcutDG{dNo}(sNo) = st.wellSol(2).wcut;
        
        st = ohDGReorder{sNo};
        statesDGReorder{dNo}{sNo} = st;
        wcutDGReorder{dNo}(sNo) = st.wellSol(2).wcut;
        
        r = rhDG{sNo};
        tDG{dNo}(sNo) = r.WallTime;
        itDG{dNo}(sNo) = cellfun(@(r) r.NonlinearReport{1}.TransportSolver.Iterations, r.StepReports);

        
        r = rhDGReorder{sNo};
        tDGReorder{dNo}(sNo) = r.WallTime;        
        itDGReorderAll{dNo}(:,sNo) = r.StepReports{1}.NonlinearReport{1}.TransportSolver.StepReports{1}.NonlinearReport{1}.Iterations;
        itDGReorder{dNo}(sNo) = mean(itDGReorderAll{dNo}(:,sNo));
        
    end
    
end

%% Plot paprams


pos = [-1000, 0, 800, 600];
azel = [115, 20];
fs = 15;
pth = fullfile(mrstPath('dg'), 'examples', 'ecmor-xvi', 'pebi', 'fig');
W(1).name = 'INJ';
W(2).name = 'PROD';
pltWell = @() plotWell(G, W, 'color', 'k', 'height', 50);

if 1
    savepng = @(name) print(fullfile(pth, name), '-dpng', '-r300');
    saveeps = @(name) print(fullfile(pth, name), '-depsc');
else
    savepng = @(name) [];
    saveeps = @(name) [];
end


%% Plot rock props

close all

% hpos = [0.2476    0.06    0.5392    0.00500]
hpos = [0.2539    0.25    0.5273    0.04];
cpos = [0.2539    0.22    0.5273    0.02];

perm = modelFI.rock.perm(:,1);
figure('position', pos)
hold on
plotCellData(G, log10(perm));
hold off; axis equal tight
ax = gca;
[ax.XTick, ax.YTick, ax.ZTick] = deal([]);
[c, h] = colorbarHist(perm, [min(perm), max(perm)], 'South', 50, true);
ax.FontSize = fs;
% colormap(jet);
h.Position = hpos;
c.Position = cpos;
view(azel)
savepng('pebi-perm');

poro = modelFI.rock.poro(:,1);
figure('position', pos)
hold on
plotCellData(G, poro);
hold off; axis equal tight
ax = gca;
[ax.XTick, ax.YTick, ax.ZTick] = deal([]);
[c, h] = colorbarHist(poro, [min(poro), max(poro)], 'South', 50, false);
ax.FontSize = fs;
% colormap(jet);
h.Position = hpos;
c.Position = cpos;
view(azel)
savepng('pebi-poro');

figure('position', pos)
hold on
plotGrid(G, 'facec', 'none');
plotGrid(G, [W.cells], 'facec', 'k');
pltWell();
hold off; axis equal off
ax = gca;
[ax.XTick, ax.YTick, ax.ZTick] = deal([]);
view(azel)
savepng('pebi-wells');

%%

close all

stateNo = [3, 10, 20];
thres = 0.6;
for dNo = 1:numel(degree)
    for sNo = 1:numel(stateNo)
        
        figure('position', pos);
        plotGrid(G, 'facec', 'none', 'edgec', [1,1,1]*0.7);
        pltWell();
        plotGrid(G, [W.cells], 'facec', 'k');
        s = statesDGReorder{dNo}{stateNo(sNo)}.s(:,1);
        keep = s > thres;
        plotCellData(G, s(keep), keep);
        hold off; axis equal off
        ax = gca;
        [ax.XTick, ax.YTick, ax.ZTick] = deal([]);
        view(azel)
        caxis([thres,1])
        savepng(['pebi-dg', num2str(degree(dNo)), '-sat', num2str(sNo)]);
        
    end
end

%%

close all

pp =  [-1000, 0, 500, 800];
fss = 20;

figure('Position', pp);
colorbar();
axis off
caxis([thres, 1]);
ax = gca;
ax.FontSize = fss;
savepng('pebi-satbar')

%%

close all
clr = lines(2);
dgarg = {'linewidth', 6};
dgrarg = {'linewidth', 2};
dtt = cumsum(dtvec)/day;
figure
hold on
for dNo = 1:numel(degree)
    plot(dtt, wcutDG{dNo}, '--', dgarg{:}, 'color', clr(dNo,:))
    plot(dtt, wcutDGReorder{dNo}, '-', dgrarg{:}, 'color', clr(dNo,:))
    axis([0 dtt(end) 0 1])
    box on
end
legend({'dG(0)', 'dG(0), reordered', 'dG(1)', 'dG(1), reordered'}, 'location', 'northwest')
hold off
ax = gca;
ax.FontSize = fs;
xlabel('Time [days]');
ylabel('Water cut');
saveeps('pebi-wcut');

%%

close all

tol = 5;
cells = abs(G.cells.centroids(:,2) - xmax(2)/2) < tol;
cells = find(cells);

cells = cells(1:5);

clrs = jet(numel(cells));
figure;
hold on

for cNo = 1:numel(cells)
    
    s0 = cellfun(@(st) st.s(cNo,1), statesDG{1});
    s1 = cellfun(@(st) st.s(cNo,1), statesDG{2});
    plot(s0, '--', 'color', clrs(cNo,:));
    plot(s1, '-', 'color', clrs(cNo,:));
    
end

%% Plot Iterations


close all
clr = lines(2);

dgclr = clr(1,:);
dgsz = 25;

dgrclr = clr(2,:);
dgrsz = 15;

dgrmclr = brighten(dgrclr, -0.5);
dgrmsz = 10;
tvec = cumsum(schedule.step.val)/year;
itmax = max(cellfun(@max, itDG));

for dNo = 1:numel(degree)
    
    figure('Position', [-1000, 0, 600, 300])
    
    hold on
%     its = [itDG{dNo}, itDGReorder{dNo}];
    plot(tvec, itDG{dNo}, '.-', 'color', dgclr, 'markerSize', dgsz)
    plot(tvec, itDGReorder{dNo}, '.-', 'color', dgrclr, 'markerSize', dgrsz)
    maxIt = max(itDGReorderAll{dNo}, [], 1);
    plot(tvec, maxIt, '.', 'color', dgrmclr, 'markerSize', dgrmsz)
    box on
    xlabel('Time [Years]')
%     ylabel('Nonlinear iterations')
    dgn = ['dG(', num2str(degree(dNo)), ')'];
    legend({dgn, [dgn, ' reordered'], [dgn, ' reordered, max']}, 'location', 'southeast')
    ylim([0, itmax+10])
    ax = gca;
    ax.YScale = 'log';
    ax.YTick = [0.01, 0.1, 1, 10, 100];
    ax.FontSize = 16;
    
    drawnow(); pause(0.1);
    
    saveeps(['pebi-dg', num2str(degree(dNo)), '-it']);
    
    
end
