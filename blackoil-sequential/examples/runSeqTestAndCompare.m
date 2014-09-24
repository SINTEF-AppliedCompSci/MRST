mrstModule add ad-unittest
% [G, rock, fluid, deck, state] = setupSPE1();
testcase = TestSimpleOW();
mrstModule add ad-fi deckformat mrst-gui ad-core ad-blackoil blackoil-sequential ad-unittest


model    = testcase.model;
state    = testcase.state0;
schedule = testcase.schedule;
rock     = testcase.rock;
G        = model.G;
state.wellSol = initWellSolAD(schedule.control(1).W, model, state);


%%


clear pressureModel
mrstModule add blackoil-sequential
pressureModel = PressureOilWaterModel(G, rock, model.fluid);
transportModel = TransportOilWaterModel(G, rock, model.fluid);



seqModel = SequentialPressureTransportModel(pressureModel, transportModel);

[ws_split, states_split, report_split] = simulateScheduleAD(state, seqModel, schedule);


%% Run the entire schedule
[ws_fi, states_fi, report_fi] = simulateScheduleAD(state, model, schedule);




%% Plot the well solutions and simulator states
% We setup interactive viewers for both well solutions and the reservoir
% states.

time = {report_fi.ReservoirTime, report_split.ReservoirTime};
ws = {ws_fi, ws_split};

plotWellSols(ws, time, 'datasetnames', {'FI', 'sequential'})

figure;
plotToolbar(G, states)
plotWell(G, schedule.control(1).W)
axis tight
view(-10, 60)