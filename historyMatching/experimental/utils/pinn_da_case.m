mrstModule add ad-blackoil ad-props ad-core mrst-gui ensemble example-suite


%https://github.com/sintefmath/pinn_deep_ritz_and_oil/blob/main/flow/hybrid_a.ipynb

G = cartGrid([50, 50], [2 2]);
G.nodes.coords = G.nodes.coords - 1;
G = computeGeometry(G);


%%

rock = makeRock(G, 1, 1);
fluid = initSimpleADIFluid('phases', 'WO', 'n', [1 1], 'mu', [1 1], 'rho', [1 1]);

model = GenericBlackOilModel(G, rock, fluid, 'gas', false);

%%
time = 0.7;
rates = [1 -1];
% check rates / well controls

% flux in/out = 1


src = [];
src = addSource(src, 1, 1, 'sat', [1 0] ); % (src, cell, val)
src = addSource(src, G.cells.num, -1, 'sat', [1 0]);

dt = repmat(time/100, 100, 1);
schedule = simpleSchedule(dt, 'src', src);

%%
state0 = initResSol(G, 1, [0 1]);

%% 

problem = packSimulationProblem(state0, model, schedule, 'PINN_DA_case');

simulatePackedProblem(problem);


%%

[wellSols, states, reports] = getPackedSimulatorOutput(problem);

%%
plotToolbar(G, states)

%%

wellInj = cellfun(@(state) state.pressure(1), states);
wellProd = cellfun(@(state) state.pressure(end), states);
figure()
hold on
plot(wellInj)
plot(wellProd)
figure()
plot(wellInj - wellProd)



