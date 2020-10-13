%% Load data

[itDG, itDGReorder, itDGReorderAll, states] = deal(cell(numel(degree),1));

for dNo = 1:numel(degree)
    
    ohDG = getOutHandler(['dg', num2str(degree(dNo))]);
    rhDG = getRepHandler(['dg', num2str(degree(dNo))]);
    
    ohDGReorder = getOutHandler(['dg', num2str(degree(dNo)), '-reorder']);
    rhDGReorder = getRepHandler(['dg', num2str(degree(dNo)), '-reorder']);

    for sNo = 1:ohDG.numelData
        
        st = ohDGReorder{sNo};
        states{dNo}{sNo} = st;
        
        rDG = rhDG{sNo};
        itDG{dNo}(sNo) = cellfun(@(r) r.NonlinearReport{1}.TransportSolver.Iterations, rDG.StepReports);
        
        rDGReorder = rhDGReorder{sNo};
        
        itDGReorderAll{dNo}(:,sNo) = rDGReorder.StepReports{1}.NonlinearReport{1}.TransportSolver.StepReports{1}.NonlinearReport{1}.Iterations;
        itDGReorder{dNo}(sNo) = mean(itDGReorderAll{dNo}(:,sNo));
    end
    
end

%% Plot paprams

fs = 17;
pos = [-1000, 0, 600, 800];
xw = G.cells.centroids([schedule.control(1).W.cells],:);
pltWell = @() plot(xw(:,1), xw(:,2), 'ok', 'markerSize', 9, 'markerFacec', 'w', 'linewidth', 2);
xwn = xw;
xwn(1,2) = xwn(1,2) + 10 ;xwn(2,2) = xwn(2,2) - 10;
pltWellNames = @() text(xwn(:,1), xwn(:,2), {'INJ', 'PROD'}, 'FontSize', fs);

pth = fullfile(mrstPath('dg'), 'examples', 'ecmor-xvi', 'spe10-channels', 'fig');

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
hpos = [0.2539    0.14    0.5273    0.04];
cpos = [0.2539    0.11    0.5273    0.02];


perm = modelFI.rock.perm(:,1);
figure('position', pos)
hold on
plotCellData(G, log10(perm), 'edgec', 'none');
pltWell(); pltWellNames();
hold off; axis equal tight
ax = gca;
ax.FontSize = fs;
[ax.XTick, ax.YTick] = deal([]);
[c, h] = colorbarHist(perm, [min(perm), max(perm)], 'South', 50, true);
% colormap(jet);
h.Position = hpos;
c.Position = cpos;
box on
savepng('spe10-perm');

poro = modelFI.rock.poro(:,1);
figure('position', pos)
hold on
plotCellData(G, poro, 'edgec', 'none');
pltWell(); pltWellNames();
hold off; axis equal tight
ax = gca;
ax.FontSize = fs;
[ax.XTick, ax.YTick] = deal([]);
[c, h] = colorbarHist(poro, [min(poro), max(poro)], 'South', 50, false);
% colormap(jet);
h.Position = hpos;
c.Position = cpos;
box on
savepng('spe10-poro');

order = states{1}{10}.order;
figure('position', pos + [0,0,-170,200])
hold on
plotCellData(G, order, 'edgec', 'none');
pltWell(); pltWellNames();
hold off; axis equal tight
ax = gca;
ax.FontSize = fs;
[ax.XTick, ax.YTick] = deal([]);
c = colorbar('location', 'southoutside');
c.Ticks = round(linspace(1, G.cells.num, 5));
% colormap(jet);
cpos = [0.15    0.13  0.73    0.02];
c.Position = cpos;
box on
savepng('spe10-order');

%% Plot saturations

close all
stateNo = [10, 20, 32];

itmax = max(max(itDGReorderAll{1}(:, stateNo)));
itmax = max(itmax, max(max(itDGReorderAll{2}(:, stateNo))));
saveNo = 0;
for dNo = 1:numel(degree)
    for sNo = stateNo

        saveNo = saveNo + 1;
        
%         figure('position', pos, 'name', ['dG(', num2str(degree(dNo)), ') step ', num2str(sNo)]);
%         hold on
%         s = states{dNo}{sNo}.s(:,1);
%         keep = s > 0.2;
%         plotGrid(G, 'facec', 'none', 'edgec', 'none');
%         plotCellData(G, s(keep), keep, 'edgec', 'none');
%         pltWell()
%         hold off; axis equal tight
%         colormap(jet)
%         caxis([0.2,0.8]);
% %         colorbar('Location', 'southoutside')
%         box on
%         ax = gca;
%         [ax.XTick, ax.YTick] = deal([]);
%         
% %         savepng(['spe10-dg', num2str(degree(dNo)), '-sat', num2str(saveNo)]);
        
        figure('position', pos, 'name', ['dG(', num2str(degree(dNo)), ') it ', num2str(sNo)]);
        hold on
        
        it = itDGReorderAll{dNo}(:,sNo);

        keep = it>0;
        plotGrid(G, 'facec', 'none', 'edgec', 'none');
        plotCellData(G, it(keep), keep, 'edgec', 'none');
        pltWell()
        hold off; axis equal tight
        colormap(jet)
        caxis([1, itmax]);
%         colorbar('Location', 'southoutside')
        box on
        ax = gca;
        [ax.XTick, ax.YTick] = deal([]);
        
        savepng(['spe10-dg', num2str(degree(dNo)), '-it', num2str(saveNo)]);

    end
end

%%

close all

pp =  [-1000, 0, 500, 800];
fss = 25;

figure('Position', pp);
colorbar();
axis off
caxis([0.2,0.8]);
colormap(jet);
ax = gca;
ax.FontSize = fss;
savepng('spe10-satbar')

figure('Position', pp);
c = colorbar();
c.Ticks = 1:itmax;
axis off
caxis([1, itmax]);
c = colormap(jet);
ax = gca;
ax.FontSize = fss;
savepng('spe10-itbar')


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
    ylabel('Nonlinear iterations')
    dgn = ['dG(', num2str(degree(dNo)), ')'];
    legend({dgn, [dgn, ' reordered'], [dgn, ' reordered, max']}, 'location', 'southeast')
    ylim([0, itmax+10])
    ax = gca;
    ax.YScale = 'log';
    ax.YTick = [0.01, 0.1, 1, 10, 100];
    ax.FontSize = 16;
    
    drawnow(); pause(0.1);
    
    saveeps(['spe10-dg', num2str(degree(dNo)), '-it']);
    
    
end


%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
