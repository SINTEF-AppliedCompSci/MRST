mrstModule add ad-core ad-blackoil ad-props
%% Set up scenario
dims = [50, 50];
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
cfl = estimateSaturationCFL(model, states{1}, max(schedule.step.val));
figure;
plotCellData(G, cfl, 'EdgeColor', 'none');
cmap = parula(100);
cmap = interp1((0.02:0.02:2), cmap, log10(2:100)');
colormap(cmap);
colorbar;

outlineCoarseGrid(G, double(cfl > 1), 'Color', 'w');
%% Set up a WENO discretization
model_weno = model;
weno = WENOUpwindDiscretization(model_weno);

%% Override the component discretization with a WENO scheme
model_weno = model_weno.validateModel();
props = model_weno.FluxDiscretization;
disp(props)
props = props.setStateFunction('FaceMobility', FaceMobility(model_weno, weno));
props = props.setStateFunction('FaceComponentMobility', FaceComponentMobility(model_weno, weno));
model_weno.FluxDiscretization = props;

[ws_weno, states_weno, report_weno] = simulateScheduleAD(state0, model_weno, schedule);
%%
model_e = model.validateModel();
fd = model_e.FluxDiscretization;
fd = fd.setFlowStateBuilder(AdaptiveImplicitFlowStateBuilder('initialStep', 0.02*day, 'verbose', true));
model_e.FluxDiscretization = fd;
[ws_e, states_e, report_e] = simulateScheduleAD(state0, model_e, schedule);

%%
model_weno_expl = model_weno;
fd = model_weno_expl.FluxDiscretization;
fd = fd.setFlowStateBuilder(AdaptiveImplicitFlowStateBuilder('initialStep', 0.02*day, 'verbose', true));
model_weno_expl.FluxDiscretization = fd;
[ws_weno_ex, states_weno_ex, report_weno_ex] = simulateScheduleAD(state0, model_weno_expl, schedule);

%% Plot saturations
G.cells.sortedCellNodes = getSortedCellNodes(G);
for i = 1:numel(states)
    figure(1); clf;
    subplot(2, 2, 1)
    plotCellData(G, states{i}.s(:, 1), 'edgecolor', 'none');
    axis equal tight
    title('FIM SPU');
    subplot(2, 2, 2)
    plotCellData(G, states_weno{i}.s(:, 1), 'edgecolor', 'none');
    axis equal tight
    title('FIM WENO')
    subplot(2, 2, 3)
    plotCellData(G, states_e{i}.s(:, 1), 'edgecolor', 'none');
    title('AIM SPU')
    axis equal tight

    subplot(2, 2, 4)
    plotCellData(G, states_weno_ex{i}.s(:, 1), 'edgecolor', 'none');
    title('AIM WENO')
    axis equal tight
end
%% Plot wells
plotWellSols({ws, ws_e, ws_weno, ws_weno_ex}, report.SimulationTime, 'datasetnames', {'FIM-SPU', 'AIM-SPU', 'FIM-WENO', 'AIM-WENO'}, 'field', 'qWs', 'SelectedWells', 2)

% <html>
% <p><font size="-1">
% Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.
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
