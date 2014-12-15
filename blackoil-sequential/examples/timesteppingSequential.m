mrstModule add ad-unittest
% [G, rock, fluid, deck, state] = setupSPE1();
% testcase = TestSPE1();
testcase = TestEGG();

mrstModule add ad-fi deckformat mrst-gui ad-core ad-blackoil blackoil-sequential ad-unittest


model    = testcase.model;
state    = testcase.state0;
schedule = testcase.schedule;
rock     = testcase.rock;
G        = model.G;
state.wellSol = initWellSolAD(schedule.control(1).W, model, state);

schedule = compressSchedule(schedule);

% lim = 20;
% schedule.step.val = schedule.step.val(1:lim);
% schedule.step.control = schedule.step.control(1:lim);
%%
timestepper = IterationCountTimeStepSelector('targetIterationCount', 5,...
                                         'minRelativeAdjustment', sqrt(eps),...
                                         'maxRelativeAdjustment', inf, ...
                                         'firstRampupStep',       1*day, ...
                                         'verbose', true);
% Instansiate a non-linear solver with the timestep class as a
% construction argument.
solver = NonLinearSolver('timeStepSelector', timestepper, ...
                            'maxiterations', 25);


clear pressureModel transportModel seqModel
mrstModule add blackoil-sequential

if isa(model, 'TwoPhaseOilWaterModel')
    pressureModel  = PressureOilWaterModel(G, rock, model.fluid);
    transportModel = TransportOilWaterModel(G, rock, model.fluid, 'nonlinearTolerance', 1e-8);
else
    pressureModel  = PressureBlackOilModel(G, rock, model.fluid);
    transportModel = TransportBlackOilModel(G, rock, model.fluid, 'nonlinearTolerance', 1e-8, ...
        'conserveWater', false, 'conserveOil', true, 'conserveGas', true);
end

mrstVerbose on
seqModel = SequentialPressureTransportModel(pressureModel, transportModel);

timer = tic();
[ws_split, states_split, report_split] = simulateScheduleAD(state, seqModel, schedule,...
                                        'NonLinearSolver', solver, 'OutputMinisteps', true);
t_split = toc(timer);


%% Run the entire schedule
solver.timeStepSelector.reset();

timer = tic();
[ws_fi, states_fi, report_fi] = simulateScheduleAD(state, model, schedule,...
                        'NonLinearSolver', solver, 'OutputMinisteps', true);
t_fi = toc(timer);



%% Plot the well solutions and simulator states
% We setup interactive viewers for both well solutions and the reservoir
% states.

time = {report_fi.ReservoirTime, report_split.ReservoirTime}; 
ws = {ws_fi, ws_split};

plotWellSols(ws, time, 'datasetnames', {'FI', 'sequential'})
