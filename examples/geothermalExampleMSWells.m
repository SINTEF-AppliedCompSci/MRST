mrstModule add test-suite
mrstModule add ad-core ad-props
mrstModule add geothermal compositional
mrstModule add multiphysics-model wellpaths upr
mrstModule add mrst-gui
mrstModule add ecmor-2022-geothermal
mrstVerbose on

%%
close all
cartDims = [100,1,100];
test    = TestCase('fractured_3d_slice_geothermal', 'useWellboreModel', false, 'cartDims', cartDims, 'dt', 0.5*hour);
test = convertToWBMultiModel(test);
testRef = TestCase('fractured_3d_slice_geothermal', 'useWellboreModel', false, 'cartDims', cartDims, 'dt', 0.5*hour);

%%
lsol = MultiPhysicsLinearSolver(test.model);
lsol.solveSubproblems = true;
nls = NonLinearSolver();
nls.useLinesearch = true;
nls.LinearSolver = lsol;
% test.model.submodels.Wellbore.parentModel.operators.T = test.model.submodels.Wellbore.parentModel.operators.T*100;
problem = test.getPackedSimulationProblem('LinearSolver', lsol);
simulatePackedProblem(problem, 'restartStep', 1);

%%
problemRef = testRef.getPackedSimulationProblem();
simulatePackedProblem(problemRef, 'restartStep', 1);

%%
close all
[~   , states   , reports   ] = getPackedSimulatorOutput(problem);
wellSols = getWellSolsFromWBState(test.model, states);
[wellSolsRef, statesRef, reportsRef] = getPackedSimulatorOutput(problemRef);

test.plot(states);
test.plot(statesRef);
dstates = cellfun(@(st1, st2) compareStructs(st1.Reservoir, st2, 'includeStructs', false, 'relative', true), states, statesRef, 'UniformOutput', true);
test.plot(dstates);
plotWellSols({wellSols, wellSolsRef});

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
