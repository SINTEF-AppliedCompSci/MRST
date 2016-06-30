%%
mrstModule add ad-core ad-blackoil blackoil-sequential spe10
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

%%
fluid = initSimpleADIFluid('phases', 'WG', 'c', [1/barsa, 2/barsa], 'cR', 2),

%% Compute solution with refined time steps
% We will now compute a solution with refined time steps. As the time-steps
% become smaller, the solution becomes more accurate. In order to achieve
% increased accuracy without manually changing the timesteps, we can use
% a automatic timestep selector based on saturation change targets. We let
% the solver aim for a maximum saturation change of 1% in each cell during
% the timesteps to get very small steps.
stepSel = StateChangeTimeStepSelector(...
          'targetProps', 's',...            % Saturation as change target
          'targetChangeAbs', 0.01,...       % Target change of 0.01
          'firstRampupStepRelative', 0.01); % Initial rampup step is dt0/100
solver = NonLinearSolver('timeStepSelector', stepSel);
% Simulate with small timesteps. As the resulting timesteps will be very
% small, it will take some time (about six minutes on the workstation where
% the example was written).
[wsFine, statesFine, repFine] = simulateScheduleAD(state, model, schedule, ...
                        'nonlinearsolver', solver, 'outputMinisteps', true);
