mrstModule add ad-unittest
% [G, rock, fluid, deck, state] = setupSPE1();
testcase = TestSPE1();
% testcase = TestEGG();

mrstModule add ad-fi deckformat mrst-gui ad-core ad-blackoil blackoil-sequential ad-unittest


model    = testcase.model;
state    = testcase.state0;
schedule = testcase.schedule;
rock     = testcase.rock;
G        = model.G;


% lim = 10;
% schedule.step.val = schedule.step.val(1:lim);
% schedule.step.control = schedule.step.control(1:lim);

%%
mrstModule add agmg
solver = NonLinearSolver('enforceResidualDecrease', false, 'useRelaxation', true);

clear pressureModel transportModel seqModel
mrstModule add blackoil-sequential


amgSolver = AGMGSolverAD();
mrstVerbose on
seqModel = getSequentialModelFromFI(model, 'pressureLinearSolver', amgSolver);

model.extraWellSolOutput = true;
seqModel.pressureModel.extraWellSolOutput = true;
seqModel.transportModel.extraWellSolOutput = true;

seqModel.stepFunctionIsLinear = true;

timer = tic();
[ws_split, states_split, report_split] = simulateScheduleAD(state, seqModel, schedule, 'NonLinearSolver', solver);
t_split = toc(timer);


%% Run the entire schedule
cprsolver = CPRSolverAD('ellipticSolver', amgSolver);

timer = tic();
[ws_fi, states_fi, report_fi] = simulateScheduleAD(state, model, schedule, 'LinearSolver', cprsolver);
t_fi = toc(timer);

%% Run schedule with outer loop enabled
seqModel.stepFunctionIsLinear = false;
seqModel.outerTolerance = 1e-5;
seqModel.outerCheckWellConvergence = false;
timer = tic();
[ws_outer, states_outer, report_outer] = simulateScheduleAD(state, seqModel, schedule, 'NonLinearSolver', solver);
t_outer = toc(timer);


%% Plot the well solutions and simulator states
% We setup interactive viewers for both well solutions and the reservoir
% states.

time = {report_fi.ReservoirTime, report_split.ReservoirTime, report_outer.ReservoirTime}; 
ws = {ws_fi, ws_split, ws_outer};

plotWellSols(ws, time, 'datasetnames', {'FI', 'sequential', 'outerloop'})
%%
for i = 1:2
    if i == 1
        states = states_fi;
        t = 'fully implicit';
    else
        states = states_split;
        t = 'sequential';
    end
    
    figure;
    plotToolbar(G, states)
    plotWell(G, schedule.control(1).W)
    axis tight
    view(-10, 60)
    title(t)
end

%%
figure; plotToolbar(G, states_fi{1}.s - states_split{1}.s)

%%
for i = 1:numel(states_split)
    states_split{i}.v = faceFlux2cellVelocity(G, states_split{i}.flux(:, 1));
end
figure; plotToolbar(G, states_split); axis tight off; colorbar
