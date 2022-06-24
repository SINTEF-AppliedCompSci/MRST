mrstModule add test-suite
mrstModule add ad-core ad-props
mrstModule add geothermal compositional
mrstModule add multiphysics-model wellpaths upr
mrstModule add mrst-gui
mrstVerbose on

%%
close all
test = TestCase('fractured_3d_slice_geothermal', 'useWellboreModel', true);
test.plot(test.state0);

%%
% test.model.submodels.Reservoir.applyResidualScaling = true;
% test.model.submodels.Wellbore.parentModel.applyResidualScaling = true;
lsol = MultiPhysicsLinearSolver(test.model);
problem = test.getPackedSimulationProblem('LinearSolver', lsol);
simulatePackedProblem(problem);

%%
[wellSols, states, reports] = getPackedSimulatorOutput(problem);

test.plot(states);
plotWellSols(wellSols);

%%
traj = linspace(0,max(test.model.G.nodes.coords(:,3)), 100)';
traj = [repmat(min(test.model.G.cells.centroids(:,1:2), [], 1), 100, 1), traj];
% traj = flipud(traj);

W = MultisegmentWellModel(test.model, traj);

wellState = initResSol(W.parentModel.G, 1, 1);
W = W.validateModel;
wellState = W.validateState(wellState);

schedule = simpleSchedule(rampupTimesteps(10*day, 1*day));

[~, states] = simulateScheduleAD(wellState, W, schedule);

%%

mm = MultiPhysicsModel({test.model, W});

models = mm.getModels();
names  = mm.getModelNames();
state0 = struct();
schedule = simpleSchedule(rampupTimesteps(10*day, 1*day));
for i = 1:mm.numModels()
    model = models{i};
    name = names{i};
    state0.(name) = initResSol(model.G, 1, 1);
    schedule.control(1).(names{i}) = schedule.control(1);
end
schedule.control = rmfield(schedule.control, {'W', 'src', 'bc'});
simulateScheduleAD(state0, mm, schedule);
