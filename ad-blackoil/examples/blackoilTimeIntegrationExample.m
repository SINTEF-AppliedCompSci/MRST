mrstModule add ad-core ad-blackoil deckformat ad-props linearsolvers
mrstVerbose on
[G, rock, fluid, deck, state] = setupSPE1();
[state0, model, schedule] = initEclipseProblemAD(deck);

%% Fully-implicit
% The default discretization in MRST is fully-implicit. Consequently, we
% can use the model as-is.
implicit = packSimulationProblem(state0, model, schedule, 'SPE1_ex', 'Name', 'Fully-Implicit');
%% Explicit solver
% This solver has a time-step restriction based on the CFL condition in
% each cell. The solver estimates the time-step before each solution.
model_explicit = model.validateModel();
% Get the flux discretization
flux = model_explicit.FluxDiscretization;
% Use explicit form of flow state
fb = ExplicitFlowStateBuilder();
flux = flux.setFlowStateBuilder(fb);
model_explicit.FluxDiscretization = flux;
explicit = packSimulationProblem(state0, model_explicit, schedule, 'SPE1_ex', 'Name', 'Explicit');

%% Adaptive implicit
% We can solve some cells implicitly and some cells explicitly based on the
% local CFL conditions. For many cases, this amounts to an explicit
% treatment far away from wells or driving forces. The values for
% estimated composition CFL and saturation CFL to trigger a switch to
% implicit status can be adjusted.
model_aim = model.validateModel();
flux = model_aim.FluxDiscretization;
fb = AdaptiveImplicitFlowStateBuilder();
flux = flux.setFlowStateBuilder(fb);
model_aim.FluxDiscretization = flux;
aim = packSimulationProblem(state0, model_aim, schedule, 'SPE1_ex', 'Name', 'Adaptive-Implicit (AIM)');
%% Simulate the problems
problems = {implicit, explicit, aim};
simulatePackedProblem(problems);
%% Get output and plot the well results
% There are oscillations in the well curves. Increasing the CFL limit
% beyond unity will eventually lead to oscillations.
[ws, states, reports, names, T] = getMultiplePackedSimulatorOutputs(problems);
plotWellSols(ws, T, 'datasetnames', names);
%% Plot the fraction of cells which are implicit
imp_frac = cellfun(@(x) sum(x.implicit)/numel(x.implicit), states{3});
figure;
plot(100*imp_frac);
title('Number of implicit cells in AIM')
ylabel('Percentage of implicit cells')
xlabel('Step index')
%% Plot the time-steps taken
% The implicit solvers use the control steps while the explicit solver must
% solve many substeps.
figure; hold on
markers = {'o', '.', 'x'};
for i = 1:numel(reports)
    dt = getReportMinisteps(reports{i});
    x = (1:numel(dt))/numel(dt);
    plot(x, dt/day, 'marker', markers{i})
end
xlabel('Timestep length [day]')
legend(names)