mrstModule add ad-core ad-blackoil ad-props
%% Set up scenario
dims = [15, 15];
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

n = 100;
dt = repmat(time/n, n, 1);
schedule = simpleSchedule(dt, 'W', W);

state0 = initResSol(G, 100*barsa, [0, 1]);
%% Simulate base case
[ws, states, report] = simulateScheduleAD(state0, model, schedule);
%% Set up a WENO discretization
model_weno = model;
weno = WENOUpwindDiscretization(model_weno);

%% Override the component discretization with a WENO scheme
model_weno = model_weno.validateModel();
props = model_weno.FluxDiscretization;
disp(props)
props.FaceComponentMobility = FaceComponentMobility(model_weno, weno);
props.FaceMobility = FaceMobility(model_weno, weno);

model_weno.FluxDiscretization = props;

[ws_weno, states_weno, report_weno] = simulateScheduleAD(state0, model_weno, schedule);
%% Plot saturations
for i = 1:numel(states)
    figure(1); clf;
    subplot(1, 2, 1)
    plotCellData(G, states{i}.s(:, 1));
    title('SPU');
    subplot(1, 2, 2)
    plotCellData(G, states_weno{i}.s(:, 1));
    title('WENO')
end
%% Plot wells
plotWellSols({ws, ws_weno}, report.SimulationTime, 'datasetnames', {'SPU', 'WENO'}, 'field', 'qWs', 'SelectedWells', 2)