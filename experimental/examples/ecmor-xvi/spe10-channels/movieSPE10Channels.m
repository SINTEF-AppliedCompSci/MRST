mrstModule add project-postprocess

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
        states{dNo}{sNo}.iterations = rDGReorder.StepReports{1}.NonlinearReport{1}.TransportSolver.StepReports{1}.NonlinearReport{1}.Iterations;
        itDGReorder{dNo}(sNo) = mean(itDGReorderAll{dNo}(:,sNo));
    end
    
end

%%

xw = G.cells.centroids([schedule.control(1).W.cells],:);
pltWell = @() plot(xw(:,1), xw(:,2), 'ok', 'markerSize', 9, 'markerFacec', 'w', 'linewidth', 2);
xwn = xw;
xwn(1,2) = xwn(1,2) + 10 ;xwn(2,2) = xwn(2,2) - 10;
pltWellNames = @() text(xwn(:,1), xwn(:,2), {'INJ', 'PROD'}, 'FontSize', fs);

%%

close all

fig = figure('position', [-1000, 0, 1000, 600]);

itMax = max(cellfun(@(s) max(s.iterations), states{1}));
itMax = max(itMax, max(cellfun(@(s) max(s.iterations), states{2})));

cmapIt = [1,1,1; jet];


ca = [0.2, 0.8];
azel = [0,90];

subplot(1,4,1)
hold on
h{1} = plotCellData(G, states{1}{1}.s(:,1), 'edgec', 'none');
caxis(ca);
view(azel);
axis equal tight
ax = gca;
[ax.XTickLabel, ax.YTickLabel] = deal([]);
box on
colormap(jet)
pltWell();
hold off

subplot(1,4,2)
hold on
h{2} = plotCellData(G, states{1}{1}.iterations, 'edgec', 'none');
caxis([0,itMax]);
view(azel);
axis equal tight
ax = gca;
[ax.XTickLabel, ax.YTickLabel] = deal([]);
box on
colormap(cmapIt)
pltWell();
hold off

subplot(1,4,3)
hold on
h{3} = plotCellData(G, states{2}{1}.s(:,1), 'edgec', 'none');
caxis(ca);
view(azel);
axis equal tight
ax = gca;
[ax.XTickLabel, ax.YTickLabel] = deal([]);
box on
pltWell();
hold off

subplot(1,4,4)
hold on
h{4} = plotCellData(G, states{2}{1}.iterations, 'edgec', 'none');
caxis([0,itMax]);
view(azel);
axis equal tight
ax = gca;
[ax.XTickLabel, ax.YTickLabel] = deal([]);
box on
pltWell();
hold off

stt = {states{1}, states{1}, states{2}, states{2}};

location = fullfile(mrstPath('dg'), 'examples', 'ecmor-xvi', 'spe10-channels', 'fig');

makeMovieFromState(G, stt, {'s', 'iterations', 's', 'iterations'}, [1,1,1,1], fig, h, 'location', location, 'duration', 15);


%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
