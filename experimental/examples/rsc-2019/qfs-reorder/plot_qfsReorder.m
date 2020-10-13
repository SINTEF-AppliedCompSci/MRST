%% Load results

[ws, states, reports] = getMultiplePackedSimulatorOutputs(problems);

%% Common parameters for plotting

% Plot wells
pw = @(G,W) plot(G.cells.centroids([W.cells], 1)     , ...
                 G.cells.centroids([W.cells], 2)     , ...
                 'ok', 'markerSize', 8, 'markerFaceColor', 'w', 'lineWidth', 2);

% Figures 
pos  = [-1000, 0, 500, 500];
posv = [-1000, 0, 500, 500];
fontSize = 12;
pth = fullfile(mrstPath('dg'), 'examples', 'rsc-2019', 'qfs-reorder', 'fig');
if 1
    savepng = @(name) print(fullfile(pth, name), '-dpng', '-r300');
    saveeps = @(name) print(fullfile(pth, name), '-depsc');
else
    savepng = @(name) [];
    saveeps = @(name) [];
end

hpos = [0.1300 0.1146 0.7750 0.0727];
cpos = [0.1300 0.07 0.7750 0.03];

% colors
gr = [1,1,1]*0.8;
cmap = parula;
frac = 0.5;
cmap = cmap*frac + (1-frac);

%% Plot saturation

names = {'reorder', 'gravity'};

% close all
frac = 0.8;

refFactor = [];
for mNo = 2
    for tNo = timeSteps
        figure('Position', pos);
        hold on
        st = states{mNo}{tNo}; 
        plotCellData(G, states{mNo}{tNo}.order);
        c = sum(order == order',2) > 1;
        plotGrid(G, c, 'facec', 'white');
        pw(G, W)
        s = states{mNo}{tNo}.s(:,1);
        unstructuredContour(G, s, 7, 'color', 'k', 'linew', 2);
        order = st.order;
        hold off
        axis equal tight; box on;
        ax = gca;
        [ax.XTickLabel, ax.YTickLabel] = deal({});
        colormap(jet)
%         savepng(['qfs-', names{mNo}, '-', num2str(tNo)]);

    end
end

%% Spyplots

close all
timeSteps = [5, 15, 30, 50, 85];
for mNo = 2
    setup = problems{mNo}.SimulatorSetup;
    model = setup.model.transportModel;
    drivingForces = setup.schedule.control;
    for tNo = timeSteps

        nc = G.cells.num;
        st  = states{mNo}{tNo}; 
        st0 = states{mNo}{tNo-1}; 
        
        st = model.initStateFunctionContainers(st);
        p = model.getEquations(st0, st, dt, drivingForces, ...
                      'iteration', 0, 'resOnly', false);
        A = p.getLinearSystem();
        
        
%         o0 = getTopologicalFluxOrdering(G, st, true);
        o0 = st.order;
        [os0, ixSort] = sort(o0);
        o = (1:G.cells.num);
        o = o(ixSort);
      
        P = sparse((1:nc), o, 1, nc, nc);
        Ar = P*A*P';
        figure('Position', pos, 'name', [names{mNo}, ' reordered']);
        [ii, jj] = find(Ar);
        plot(jj, ii, '.k', 'markerSize', 12);
        
        [os,n] = rlencode(os0);
        cPos = [0;cumsum(n)]+1;
        
        hold on
        cycles = os(n>1);
        clr = lines(numel(cycles));
        for cNo = 1:numel(cycles)
            cycle = cycles(cNo);
            ix = (cPos(cycle):cPos(cycle+1)-1)';
            x = [ix(1)-0.5, ix(1)-0.5;
                 ix(1)-0.5, ix(end)+0.5;
                 ix(end)+0.5, ix(end)+0.5;
                 ix(end)+0.5, ix(1)-0.5;
                 ix(1)-0.5, ix(1)-0.5];
            plot(x(:,1), x(:,2), 'color', clr(cNo,:), 'linew', 1.5);
            
            ind = ii >= ix(1) & ii <= ix(end) & jj >= ix(1) & jj <= ix(end);
            plot(jj(ind), ii(ind), '.', 'markerSize', 12, 'color', clr(cNo,:));
        end
        hold off
        axis equal
        axis([-0.5, nc+0.5, -0.5, nc+0.5]); 
        ax = gca;
        [ax.XTick, ax.YTick] = deal([]);
        ax.YDir = 'reverse';
        drawnow(); pause(0.5);
        saveeps(['qfs-spy-re-', names{mNo}, '-', num2str(tNo)])
        
        
        figure('Position', pos, 'name', names{mNo});
        [ii, jj] = find(A);
        plot(jj, ii, '.k', 'markerSize', 12);
    
        hold on
        for cNo = 1:numel(cycles)
            cycle = cycles(cNo);
            ix = (cPos(cycle):cPos(cycle+1)-1)';
            ind = any(ii == o(ix),2) & any(jj == o(ix),2);
            plot(jj(ind), ii(ind), '.', 'markerSize', 12, 'color', clr(cNo,:));
        end
        hold off
        
        axis equal
        axis([-0.5, nc+0.5, -0.5, nc+0.5]); 
        ax = gca;
        [ax.XTick, ax.YTick] = deal([]);
        ax.YDir = 'reverse';
        drawnow(); pause(0.5);
        saveeps(['qfs-spy-', names{mNo}, '-', num2str(tNo)])
        
        figure('Position', pos);
        hold on
        st = states{mNo}{tNo}; 
        plotGrid(G, 'facec', 'none');
        for cNo = 1:numel(cycles)
            cycle = cycles(cNo);
            ix = cPos(cycle):cPos(cycle+1)-1;
            c = any(st.order == os0(ix)',2);
%             c = any(st.order == ixSort(ix)',2);
            plotGrid(G, c, 'facec', clr(cNo,:));
        end
        pw(G, W)
        s = states{mNo}{tNo}.s(:,1);
        unstructuredContour(G, s, 7, 'color', 'k', 'linew', 2);
        order = st.order;
        hold off
        axis equal tight; box on;
        ax = gca;
        [ax.XTickLabel, ax.YTickLabel] = deal({});
        colormap(cmap)
        drawnow(); pause(0.5);
        saveeps(['qfs-sat-', names{mNo}, '-', num2str(tNo)])
%         for cNo = 1:nc
%             x = G.cells.centroids(cNo,:);
%             text(x(:,1), x(:,2), num2str(cNo));
%         end
%         savepng(['qfs-', names{mNo}, '-', num2str(tNo)]);
        

    end
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
