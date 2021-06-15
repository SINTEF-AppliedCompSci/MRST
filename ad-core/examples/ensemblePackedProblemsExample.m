%% Managing multiple realizations as packed problems
% We generate a 1D Buckley-Leverett displacement with 50 different porosity
% realizations. The permeability is correlated with the porosity.
mrstModule add ad-core ad-blackoil ad-props mrst-gui
rng(0);
n = 50;
G = cartGrid(100, 1000*meter);
G = computeGeometry(G);

% One pore-volume with maximum porosity
time = 5*year;
irate = sum(G.cells.volumes)*0.5/time;

poro = zeros(G.cells.num, n);
for i = 1:n
    poro(:, i) = reshape(gaussianField(G.cartDims, [0.01 0.5], [11 3], 2.5), [], 1);
end
K = poro.^3.*(1e-5)^2./(0.81*72*(1-poro).^2);
bc = [];
bc = fluxside(bc, G, 'XMin', irate, 'sat', [1, 0]);
bc = pside(bc, G, 'XMax', 100*barsa, 'sat', [1, 0]);

dt = rampupTimesteps(time, time/100);
schedule = simpleSchedule(dt, 'bc', bc);
fluid = initSimpleADIFluid('mu', [1, 5]*centi*poise, 'n', [2, 2], 'phases', 'wo', 'cR', 1e-8/barsa, 'c', [1e-5, 1e-4]/barsa);

state0 = initResSol(G, 100*barsa, [0, 1]);

baseName = 'EnsembleExample1D';
%% Plot a few realizations
% Plot the porosity and permeability of a few realizations
subs = 1:min(4, n);
x = G.cells.centroids;
l = arrayfun(@(x) sprintf('Realization %d', x), subs, 'UniformOutput', false);
clf;
yyaxis left
plot(x, poro(:, subs));
ylim([0, 0.6])
legend(l, 'location', 'north', 'Orientation', 'horizontal')
ylabel('Porosity')
xlabel('Position')
yyaxis right
plot(x, K(:, subs)/(milli*darcy));
legend(l, 'location', 'north', 'Orientation', 'horizontal')
ylabel('Permeability [md]')
xlabel('Position')

%% Build a model and corresponding problem for each realizations
problems = cell(1, n);
for poroNo = 1:n
    k = K(:, poroNo);
    p = poro(:, poroNo);
    caseName = sprintf('Case %d', poroNo);
    rock = makeRock(G, k, p);
    model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    description = sprintf('Average porosity %1.2f, average perm %1.2f md', mean(p), mean(k)/(milli*darcy));
    problems{poroNo} = packSimulationProblem(state0, model, schedule, baseName, 'Name', caseName, 'Description', description);
end
% Set up a problem manager
ppm = PackedProblemManager(problems);

%% Simulate in background with maximum number of threads
ppm.simulateProblemsBatch();
%% Get the output
% There is obviously a lot of output. We could set the 'readFromDisk' flag
% to false to avoid reading all the reservoir results into memory.
[ws, states, reports, names] = getMultiplePackedSimulatorOutputs(problems);
%% Plot the saturation distribution for the entire ensemble 
% The color corresponds to the realization number
nstep = numel(schedule.step.val);
np = numel(states);
colors = parula(np);

figure(1); clf;
colormap(colors);
for stepNo = 1:nstep
    clf; hold on;
    avg = 0;
    for i = 1:np
        s = states{i}{stepNo}.s(:, 1);
        plot(x, s, 'Color', [colors(i, :), 0.4], 'linewidth', 2);
        avg = avg + s;
    end
    plot(x, avg/np, 'k', 'linewidth', 2);
    colorbar
    ylabel('S_w')
    caxis([1, np])
    drawnow
end
%% Plot the water breakthrough for the ensemble
figure(1); clf;
clf; hold on;
for i = 1:np
    v = cellfun(@(x) x.s(end, 1), states{i});
    ph = plot(v, 'Color', [colors(i, :), 0.4], 'linewidth', 2);
end
colorbar
caxis([1, np])
title('Water breakthrough')

%% Plot the distribution of the front as a function of time
tol = 0.01;
bins = linspace(1, G.cells.num+1, 25);
for stepNo = 1:nstep
    clf; hold on;
    positions = zeros(np, 1);
    mv = G.cells.num;
    Mv = 0;
    for i = 1:np
        s = states{i}{stepNo}.s(:, 1);
        pos = find(s >= tol, 1, 'last');
        if ~isempty(pos)
            positions(i) = pos;
            Mv = max(Mv, pos);
            mv = min(mv, pos);
        end
    end
    histogram(positions - 0.5, bins);
    ylim([0, 10]);
    title(sprintf('Step %d: Front between %d and %d', stepNo, mv, Mv));
    drawnow
    pause(0.5);
end

%% Copyright Notice
%
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
