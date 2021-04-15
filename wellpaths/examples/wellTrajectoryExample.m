%% Define a model
mrstModule add ad-core ad-blackoil diagnostics wellpaths
km = kilo*meter;
pdims = [5, 4, 0.01]'*km;

% dims = [50, 40, 10];
dims = [25, 20, 5];

G = cartGrid(dims, pdims);
G = computeGeometry(G);
%% Define a forked well
% We define four curves, with some common points to combine them into a
% single well path.

p1 = [0.25, 0.5, 0.9, 1.5]'*km;

n = numel(p1);

z = [.25, 0.5, 0.75, .80]'*pdims(3);

a = [p1(1)*ones(n, 1), p1];
b = [p1,               p1(1)*ones(n, 1)];
c = [p1,               p1];

origin = [0, 0, 0; a(1, :), z(1)];

% First foru well paths.
wellpath0 = makeSingleWellpath(origin);
wellpath1 = makeSingleWellpath([a, z]);
wellpath2 = makeSingleWellpath([b, z]);
wellpath3 = makeSingleWellpath([c, z]);
% Combine into single well path
wellpath_fork = combineWellPaths({wellpath0, wellpath1, wellpath2, wellpath3});

% Plotting
clf;
plotWellPath(wellpath_fork);
%% Make a comb-like well
x = [3, 3.5, 4, 4.5]'*km;
y = [3.5, 3, 2.5, 2]'*km;

y0 = [3.75*km; y];
x0 = x(end)*ones(numel(y)+1, 1);

dz = 1/numel(y);
z0 = (0:dz:1)'*pdims(3);

wp0 = [x0, y0, z0];

n0 = numel(x0);

wellpaths = cell(5, 1);

wellpaths{1} = makeSingleWellpath(wp0);

for i = 1:4
    n = numel(x);
    wp = [x, repmat(y(i), n, 1), repmat(z0(i+1), n, 1)];
    wellpaths{i+1} = makeSingleWellpath(wp(end:-1:1, :));
end

wellpath_comb = combineWellPaths(wellpaths);

%% Determine the cells
[cells_fork, segInd_fork, ~, ~, DT] = findWellPathCells(G, wellpath_fork);
[cells_comb, segInd_comb] = findWellPathCells(G, wellpath_comb, 'triangulation', DT);

%% Plot the well trajectories
clf;
plotWellPath(wellpath_comb);
plotWellPath(wellpath_fork);

plotGrid(G, 'facec', 'none', 'edgea', .1)
view(40, 56)
%% Plot wells individually + realized cell perforations
close all
for i = 1:2
    if i == 1
        wp = wellpath_comb;
        si = segInd_comb;
        c = cells_comb;
        v = [30, 12];
    else
        wp = wellpath_fork;
        si = segInd_fork;
        c = cells_fork;
        v = [160, 40];
    end
    
    figure;
    subplot(2, 1, 1)
    plotWellPath(wp);
    view(v);
    title('Well trajectory with segment indicators')
    axis tight off

    subplot(2, 1, 2)
    plotCellData(G, si, c)
    view(v);
    title('Cells with segment indicators')
    axis tight off
end
%% Set up actual simulation wells from the well paths
% Initial reservoir state
initSat = [.1, .9];
state = initResSol(G, 200*barsa, initSat);
% Rock
rock = makeRock(G, 300*milli*darcy, 0.5);

time = 10*year;
rate = .25*sum(poreVolume(G, rock))/time;

segInd = cell(2, 1);
W = [];
[W, segInd{1}]  = getWellFromPath(W, G, rock, wellpath_fork, ...
    'comp_i', [1 0], 'val', 300*barsa, 'type', 'bhp', 'sign', -1, 'Name', 'Producer');
[W, segInd{2}] = getWellFromPath(W, G, rock, wellpath_comb,...
    'comp_i', [1 0], 'val', rate, 'type', 'rate', 'sign', 1, 'Name', 'Injector');

%% Set up simulation model
mrstModule add ad-core ad-blackoil ad-props

fluid = initSimpleADIFluid('rho', [1000, 700, 100]*kilogram/meter^3, ...
                           'mu', [1, 10 1]*centi*poise);

model = TwoPhaseOilWaterModel(G, rock, fluid);

model.extraStateOutput = true;
model.extraWellSolOutput = true;
%% Build a schedule
n = 50;
dt = time/n;
timesteps = repmat(dt, n, 1);
schedule = simpleSchedule(timesteps, 'W', W);

[ws, states] = simulateScheduleAD(state, model, schedule);
%% Plot well curves + reservoir properties
mrstModule add mrst-gui
figure;
plotToolbar(G, states);
axis tight off
view(40, 56);
T = cumsum(schedule.step.val);

plotWellSols(ws, T)
%% Plot the saturation front
close all
plotWellPath(wellpath_comb);
plotWellPath(wellpath_fork);

plotGrid(G, 'facec', 'none', 'edgea', .2);
view(-60, 60);
h = nan;
for i = 1:numel(states)
    if ishandle(h); delete(h); end
    
    s = states{i}.s(:, 1);
    h = plotCellData(G, s, s > 0.2, 'edgea', .2);
    title(formatTimeRange(T(i)));
    pause(0.25)
end
%% Apply some flow diagnostics...
diagstate = states{end};
[state_split, Wc] = expandWellCompletions(diagstate, W, segInd);
D  = computeTOFandTracer(state_split, G, rock, 'wells', Wc);
WP = computeWellPairs(state_split, G, rock, Wc, D);
%% Plot results
figure; plotToolbar(G, D)
plotWellPath(wellpath_comb);
plotWellPath(wellpath_fork);

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
