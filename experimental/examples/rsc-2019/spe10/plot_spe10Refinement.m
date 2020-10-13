%% Load results

[ws, states, r] = deal(cell(numel(degree), 5));

mdlIx = 1:3;
for dNo = 1:numel(degree)
    for mNo = mdlIx
        [ws{dNo, mNo}, states{dNo, mNo}, r{dNo, mNo}] ...
            = getPackedSimulatorOutput(problems{dNo, mNo}, 'readFromDisk', false);
    end
end
[wsWENO, statesWENO, reportsWENO] = getPackedSimulatorOutput(weno, 'readFromDisk', false);

%% Load Coarse 

[wsCoarse, statesCoarse, repCoarse] = deal(cell(2,1));
for dNo = 1:2
    [wsCoarse{dNo}, statesCoarse{dNo}, repCoarse{dNo}] ...
        = getPackedSimulatorOutput(coarseProblems{dNo}, 'readFromDisk', false);
end

%% Get iterations

reports = cell(size(r));
iterations = cell(numel(degree), 4);
for dNo = 1:numel(degree)
    for mNo = 1:3
        fprintf('Getting iterations for %s ... \n', names{mNo});
        n   = r{dNo,mNo}.numelData;
        rep = cell(n,1);
        for sNo = 1:n
            rep{sNo} = r{dNo,mNo}{sNo};
        end
        if contains(names{mNo}, 'reorder') && ~isempty(r{dNo, mNo})
            iterations{dNo, mNo} = getReorderingTransportIterations(rep);
        else
            for sNo = 1:numel(rep)
                fprintf('Getting iterations for state %d ... \n', sNo);
                st = states{dNo, mNo}{sNo};
                if isfield(st, 'G'); g = st.G; else; g = G; end
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

cmap = winter;
cmap = cmap(end:-1:1, :);

%% Plot Saturation profiles for dG(0), dG(1) and WENO

close all

st      = states(:, 1);
stAdapt = states(:, 3);

timeSteps = [42];
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
            
            % Plot dG profile
            figure('position', posv, 'name', ['dG(', num2str(degree(dNo)), ') adaptive']);

            unstructuredContour(G, stAdapt{dNo}{tNo}.s(:,1), 10,'linew', 2);
            hold on
            pw(G, WF);
            axis equal tight
            box on
            caxis([0.2, 0.8]);
            ax = gca;
            [ax.XTickLabel, ax.YTickLabel] = deal({});
            colormap(cmap)
            savepng(['spe10-sat-', num2str(tNo), '-dg', num2str(degree(dNo)), '-adapt']);

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

%% Plot reordering iterations

close all

rIx = find(contains(names, 'reorder')); rIx = rIx(1);
rep = reports(:, rIx);
its = iterations(:, rIx);
itMaxdG0 = max(cellfun(@(it) max(it), its{1}(timeSteps)));
itMaxdG1 = max(cellfun(@(it) max(it), its{2}(timeSteps)));
itMax = max(itMaxdG0, itMaxdG1);

timeSteps = [20, 42, 100];
for tNo = timeSteps
    for dNo = 1:2

    % Plot reordering iterations
    figure('position', posv, 'name', ['dG(', num2str(degree(dNo)), ') iterations']);

    plotGrid(G, 'facec', 'none', 'edgec', 'none')
    c = its{dNo}{tNo} > 0;
    plotCellData(G, its{dNo}{tNo}(c), c, 'edgec', 'none'); 
    hold on
    pw(G, WF);
    axis equal tight
    box on
    caxis([1, itMax]);
    ax = gca;
    [ax.XTickLabel, ax.YTickLabel] = deal({});
    colormap(jet)
    savepng(['spe10-its-', num2str(tNo), '-dg', num2str(degree(dNo))]);
    
    end
end

%% Colorbar for iterations

figure('position', [-1000,0,230,400]);
colorbar();
caxis([1,itMax]);
colormap(jet)
ax = gca;
ax.FontSize = 14;
axis off
savepng(['spe10-its-bar']);

%% Plot adaptive refinement for dG(0) and dG(1)

close all

aIx = find(contains(names, 'adapt')); aIx = aIx(1);
st  = states(:, aIx);
% rep = reports(:, aIx);
% its = iterations(:, aIx);
timeSteps = [15, 20, 30, 43];

for tNo = timeSteps
    for dNo = 1:2
        if ~isempty(st{dNo}) && st{dNo}.numelData >= tNo

            figure('position', posv, 'name', ['dG(', num2str(degree(dNo)), ')']);

%              pink = [214, 154, 153]/255;
            hold on
            plotGrid(st{dNo}{tNo}.G, 'facec', 'none', 'edgecolor', gray);
            unstructuredContour(G, st{dNo}{tNo}.s(:,1), 10, 'linewidth', 2, 'useNodeMax', true);
%             plotCellData(G, st{dNo}{tNo}.s(:,1), 'edgec', 'none'); 
            pw(G, WF);
            axis equal tight
            box on
            caxis([0.2, 0.8]);
            ax = gca;
            [ax.XTickLabel, ax.YTickLabel] = deal({});
            colormap(cmap)
            drawnow()
            pause(0.1)
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
    iimax = -Inf;
    figure('position', itPos)
    pNames = {['dG(', num2str(dNo-1), ')']           , ...
              ['dG(', num2str(dNo-1), '), reordered'], ...
              ['dG(', num2str(dNo-1), '), adaptive'] , ...
              'WENO'};
    hold on
    for mNo = 1:3
        ii = cumsum(intIts{dNo,mNo}/model.G.cells.num);
        iimax = max(max(ii), iimax);
        plot(dtt(1:numel(ii)), ii, '-', 'linew', 2)
    end
    if dNo == 2
        ii = cumsum(intItsWENO/model.G.cells.num);
        iimax = max(iimax, max(ii));
        plot(dtt, ii, '-', 'linew', 2)
    end
    hold off
    legend(pNames, 'location', 'northwest')
    box on; grid on
    axis([0, dtt(end), 0, iimax*1.1]);
    ax = gca;
    ax.FontSize = fontSize;
    xlabel('Time (days)');
    ylabel('Iterations');
    saveeps(['spe10-iterations-', num2str(dNo-1)]);
end

%% Plot well curves

[wellSols, wcut] = deal(cell(2, 4));
for dNo = 1:2
    for mNo = [1:3,5]
        if ws{dNo, mNo}.numelData > 0
            wellSols{dNo, mNo} =  {ws{dNo, mNo}{1:ws{dNo,mNo}.numelData}};
            wcut{dNo, mNo} = cellfun(@(ws) ws(2).wcut, wellSols{dNo,mNo});
        end
    end
end
wellSolsWENO = {wsWENO{1:wsWENO.numelData}};
wcutWENO     = cellfun(@(ws) ws(2).wcut, wellSolsWENO);
for dNo = 1:2
    wellSolsBase{dNo} = {wsBase{dNo}{1:wsBase{dNo}.numelData}};
    wcutBase{dNo}     = cellfun(@(ws) ws(2).wcut, wellSolsBase{dNo});
end

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
        plot(dtt(1:numel(wcut{dNo,mNo})), wcut{dNo,mNo}, mStyle{mNo}, 'linew', lw(mNo), 'markerSize', mSize);
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
        plot(dtt(1:numel(wcut{dNo,mNo})), wcut{dNo,mNo}, mStyle{mNo}, 'linew', lw(mNo), 'markerSize', mSize);
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

%% Compare states

tIx = 20;
st  = states{1,1}{tIx};
str = states{1,2}{tIx};
sd  = compareStates(st, str);


%%

mIx = 1;
dIx = 2;
st = cell(size(states));
for dNo = dIx
    for mNo = mIx
        n = states{dNo,mNo}.numelData;
        st{dNo,mNo} = cell(n,1);
        for sNo = 1:n
            fprintf('Getting state %d of %s dG(%d) ... \n', sNo, names{mNo}, degree(dNo));
            st{dNo,mNo}{sNo} = states{dNo,mNo}{sNo};
        end
    end
end


%%


close all

figure('pos', [0,0,1100, 400])

IT0 = [250 + 150/250*50;
      415/515*50;
      (515-170)/515*50];
IT1 = [600 + 105/350*200;
       70/350*200;
       225/350*200;
       400+15/350*200];
   
x = [1,2,3,5,6,7,9];
bar(x, [IT0; IT1]);

xticklabels({'Global', 'Reordered', 'Adaptive', 'Global', 'Reordered', 'Adaptive', 'WENO'})
box on
grid on
ax = gca;
ax.FontSize = 11;
saveeps('spe10-itbar')

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
