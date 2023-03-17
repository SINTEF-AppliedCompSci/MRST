mrstModule add test-suite
mrstModule add ad-core ad-props
mrstModule add geothermal compositional dfm
mrstModule add ruden-geothermal
mrstModule add composite-model wellpaths upr
mrstModule add mrst-gui
mrstModule add ecmor-2022-geothermal
mrstModule add linearsolvers
mrstVerbose on

%%
close all
test    = TestCase('fivespot_geothermal', 'cartDims', [11, 11, 12], 'dfm', true);
test.model.FacilityModel = [];
groups = [];
hot  = {'Hot'};
cold = {'Cold-1-1', 'Cold-2-1', 'Cold-2-2', 'Cold-1-2'};
groups = addFacilityGroup(groups, hot , 'name', 'HotGroup');
groups = addFacilityGroup(groups, cold, 'name', 'ColdGroup');
test = convertToWBMultiModel(test, 'groups', groups);

ctrl = test.schedule.control(1).Wellbore;
[ctrl.W.type] = deal('group');
ctrl.W(1).status = false;
ctrl.W(3).status = false;
ctrl.groups = groups;
ctrl.groups(1).type = 'rate';
ctrl.groups(1).lims = struct('bhp', 1.001*atm);
ctrl.groups(2).type = 'bhp';

val = vertcat(ctrl.W(5).val);
ctrl.groups(1).val = val;
ctrl.groups(2).val = ctrl.W(1).val;
ctrl.groups(1).T = ctrl.W(5).T;
ctrl.groups(2).T = ctrl.W(1).T;
test.schedule.control(1).Wellbore = ctrl;

ctrl = test.schedule.control(2).Wellbore;
[ctrl.W.type] = deal('group');
ctrl.groups = groups;
ctrl.groups(1).type = 'bhp';
ctrl.groups(2).type = 'rate';
val = vertcat(ctrl.W(1:4).val);
ctrl.groups(1).val = ctrl.W(5).val;
ctrl.groups(2).val = sum(val);
ctrl.groups(1).T = ctrl.W(5).T;
ctrl.groups(2).T = ctrl.W(1).T;
test.schedule.control(2).Wellbore = ctrl;

testRef = TestCase('fivespot_geothermal', 'cartDims', [11, 11, 12], 'dfm', true);

%%
lsol = MultiPhysicsLinearSolver(test.model);
% lsol = BackslashSolverAD();
% lsol = BackslashSolverAD();
% lsol.solveSubproblems = true;
nls = NonLinearSolver();
nls.useLinesearch = true;
nls.LinearSolver = lsol;
% test.model.submodels.Wellbore.parentModel.operators.T = test.model.submodels.Wellbore.parentModel.operators.T*100;
problem = test.getPackedSimulationProblem('LinearSolver', lsol);
simulatePackedProblem(problem, 'restartStep', 1);

%%
problemRef = testRef.getPackedSimulationProblem();
simulatePackedProblem(problemRef, 'restartStep', nan);

%%
close all
[~   , states   , reports   ] = getPackedSimulatorOutput(problem);
wellSols = getWellSolsFromWBState(test.model, states);
[wellSolsRef, statesRef, reportsRef] = getPackedSimulatorOutput(problemRef);

test.plot(states);
testRef.plot(statesRef);
dstates = cellfun(@(st1, st2) compareStructs(st1.Reservoir, st2, 'includeStructs', false, 'relative', true), states, statesRef, 'UniformOutput', true);
% dstates = cellfun(@(st1, st2) compareStructs(st1.Reservoir, st2, 'includeStructs', false, 'relative', false, 'fun', @(x) x), states, statesRef, 'UniformOutput', true);
test.plot(dstates);
plotWellSols({wellSols, wellSolsRef});
plotWellSols({wellSols, wellSolsRef}, test.schedule.step.val);

plotWellSols({wellSols});
plotWellSols({wellSols}, test.schedule.step.val);

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

%%
