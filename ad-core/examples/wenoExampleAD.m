mrstModule add ad-core ad-blackoil ad-props
%% Set up scenario
nx = 50;
ny = 50;
dims = [nx, ny];
pdims = [100, 100];

G = cartGrid(dims, pdims);
G = computeGeometry(G);

rock = makeRock(G, 0.1*darcy, 0.3);
fluid = initSimpleADIFluid('n', [1, 1], 'mu', [1, 1]*centi*poise, 'rho', [100, 100], 'phases', 'wo');
model = GenericBlackOilModel(G, rock, fluid, ...
                        'water', true, 'oil', true, 'gas', false);
time = 1*year;
irate = 2*sum(model.operators.pv)/time;

W = [];
W = addWell(W, G, rock, 1, 'comp_i', [1, 0], 'val', irate, 'type', 'rate');
W = addWell(W, G, rock, G.cells.num, 'comp_i', [1, 0], 'val', 50*barsa, 'type', 'bhp');

n = 300;
dt = repmat(time/n, n, 1);
schedule = simpleSchedule(dt, 'W', W);

state0 = initResSol(G, 100*barsa, [0, 1]);
%% Simulate base case
[ws, states, report] = simulateScheduleAD(state0, model, schedule);
%% Plot the CFL condition
e = 0.01;
state = states{1};
cfl = estimateSaturationCFL(model, state, max(schedule.step.val), 'forces', schedule.control(1));
figure;
plotCellData(G, log10(cfl + e), 'EdgeColor', 'none');
colormap(jet)

% cmap = jet(100);
% cmap = interp1((0.02:0.02:2), cmap, log10(2:100)');
% colormap(cmap);
cb = colorbar;
v = [e, 0.1, 1, 10, 30];
set(cb, 'Ticks', log10(v))
set(cb, 'TickLabels', v, 'FontSize', 18)

mrstModule add streamlines;
seed = (nx:nx-1:nx*ny).';
Sf = pollock(G, state, seed, 'substeps', 1);
Sb = pollock(G, state, seed, 'substeps', 1, 'reverse', true);
hf=streamline(Sf);
hb=streamline(Sb);
set([hf; hb],'Color','k');
outlineCoarseGrid(G, double(cfl > 1), 'Color', 'w');
outlineCoarseGrid(G, ones(G.cells.num, 1), 'Color', 'k');

axis equal tight off
%% Override the component discretization with a WENO scheme
model_weno = setWENODiscretization(model);
[ws_weno, states_weno, report_weno] = simulateScheduleAD(state0, model_weno, schedule);
%% Adaptive implicit SPU
model_e = setTimeDiscretization(model, 'AIM', 'initialstep', 0.02*day, 'verbose', true);
[ws_e, states_e, report_e] = simulateScheduleAD(state0, model_e, schedule);

%% Adaptive implicit WENO
model_weno_expl = setTimeDiscretization(model_weno, 'AIM', 'initialstep', 0.02*day, 'verbose', true);
[ws_weno_ex, states_weno_ex, report_weno_ex] = simulateScheduleAD(state0, model_weno_expl, schedule);

%% Plot saturations
G.cells.sortedCellNodes = getSortedCellNodes(G);
for i = 1:numel(states)
    figure(1); clf;
    subplot(2, 2, 1)
    plotCellData(G, states{i}.s(:, 1), 'edgecolor', 'none');
    title('FIM SPU');
    axis equal tight

    subplot(2, 2, 2)
    plotCellData(G, states_weno{i}.s(:, 1), 'edgecolor', 'none');
    title('FIM WENO')
    axis equal tight

    subplot(2, 2, 3)
    plotCellData(G, states_e{i}.s(:, 1), 'edgecolor', 'none');
    title('AIM SPU')
    axis equal tight
    
    subplot(2, 2, 4)
    plotCellData(G, states_weno_ex{i}.s(:, 1), 'edgecolor', 'none');
    title('AIM WENO')
    axis equal tight
end
%% Plot the fronts at a specific time
ix = 80;
x = reshape(G.cells.centroids(:, 1), G.cartDims);
y = reshape(G.cells.centroids(:, 2), G.cartDims);
x(x == min(x(:))) = min(G.nodes.coords(:, 1));
x(x == max(x(:))) = max(G.nodes.coords(:, 1));
y(y == min(y(:))) = min(G.nodes.coords(:, 2));
y(y == max(y(:))) = max(G.nodes.coords(:, 2));

for i = 1:2
    if i == 1
        impl = states;
        expl = states_e;
    else
        impl = states_weno;
        expl = states_weno_ex;
    end
    impl = impl{ix};
    expl = expl{ix};

    vi = impl.s(:, 1);
    ei = expl.s(:, 1);
    N        = 20;
    cval     = [.5*movsum(linspace(0,1,N+1),2) 1];
    figure(i); clf; hold on
    colormap(flipud([.7*winter(128).^2+.3; 1 1 1]));
    contourf(x, y, reshape(vi, G.cartDims), cval,'EdgeColor','none');
    contour(x, y, reshape(ei, G.cartDims),  cval(2:end-1), 'k')
    outlineCoarseGrid(G, ones(G.cells.num, 1), 'Color', 'k');

    axis equal tight off;
    for wno = 1:numel(W)
        c = G.cells.centroids(W(wno).cells, :);
        plot(c(1), c(2), 'kO', 'markersize', 8, 'markerFaceColor', 'r')
    end
end

%% Plot wells
plotWellSols({ws, ws_e, ws_weno, ws_weno_ex}, report.SimulationTime, 'datasetnames', {'FIM-SPU', 'AIM-SPU', 'FIM-WENO', 'AIM-WENO'}, 'field', 'qWs', 'SelectedWells', 2)

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
