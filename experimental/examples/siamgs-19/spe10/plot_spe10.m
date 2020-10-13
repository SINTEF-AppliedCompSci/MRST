%% Plot Saturation profiles for dG(0), dG(1) and WENO

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
pth = fullfile(mrstPath('dg'), 'examples', 'siamgs-19', 'fig', 'spe10');
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

%%

close all

st      = states(:, 1);
stAdapt = states(:, 3);

timeSteps = [10, 20, 30 42];
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

%%

close all
figure('position', [0,0,700,400])
get_wcut = @(ws) cellfun(@(w) w(2).wcut, ws);

wcutDG0 = get_wcut(ws{1,1});
wcutDG1 = get_wcut(ws{2,1});
wcutWENO = get_wcut(wsWENO);

dtt = cumsum(schedule.step.val)/day;
lw = 2;
clr = lines(3);
hold on
plot(dtt, wcutDG0, 'linew', lw, 'color', clr(1,:));
plot(dtt, wcutDG1, 'linew', lw, 'color', clr(2,:));
plot(dtt, wcutWENO, 'linew', lw, 'color', clr(3,:));
ylim([0,1])
box on
ax = gca;
ax.FontSize = 12;
xlabel('Time (days)')
ylabel('Water cut');
legend({'dG(0)', 'dG(1)', 'WENO'}, 'location', 'northwest')

saveeps('spe-wcut')

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
