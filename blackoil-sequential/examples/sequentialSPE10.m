%%
[state, model, schedule] = setupSPE10_AD('layers', 15);
seqModel = getSequentialModelFromFI(model);

stepSel = StateChangeTimeStepSelector('targetProps', 's', 'targetChangeAbs', 0.25, 'firstRampupStepRelative', 0.01);
solver = NonLinearSolver('timeStepSelector', stepSel);

timer = tic();
[wsSeq, statesSeq] = simulateScheduleAD(state, seqModel, schedule, 'nonlinearsolver', solver);
t_split = toc(timer);
%% Run the entire schedule
solver.timeStepSelector.reset();

timer = tic();
[wsFIMP, statesFIMP] = simulateScheduleAD(state, model, schedule, 'nonlinearsolver', solver);
t_fi = toc(timer);

%%
plotWellSols({wsSeq, wsFIMP}, cumsum(schedule.step.val), 'datasetnames', {'Sequential', 'FIMP'})