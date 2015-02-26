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


%% Load polymer fluid from example in ad-blackoil

fn    = 'ad-blackoil/examples/polymerExamples/POLYMER.DATA';
deck  = readEclipseDeck(fn);
deck  = convertDeckUnits(deck);
fluid = initDeckADIFluid(deck);

model = OilWaterPolymerModel(G, rock, fluid, 'inputdata', deck);

% Remove gass component from wells
schedule.control.W(1).compi = [1 0];
schedule.control.W(2).compi = [0 1];

% Inject polymer
schedule.control.W(1).poly = 4;

% Set initial state
state.s = state.s(:,1:2);
state.c = zeros(G.cells.num, 1);
state.cmax = zeros(G.cells.num, 1);

% Oil rel-perm from 2p OW system.
% Needed by equation implementation function 'eqsfiOWExplictWells'.
%fluid.krO = fluid.krOW;


%% Run sequential simulation

mrstModule add agmg
solver = NonLinearSolver('enforceResidualDecrease', false, ...
   'useRelaxation', true);

clear pressureModel transportModel seqModel
mrstModule add blackoil-sequential


amgSolver = AGMGSolverAD();
mrstVerbose on
seqModel = getSequentialModelFromFI(model, ...
   'pressureLinearSolver', amgSolver);

model.extraWellSolOutput = true;
seqModel.pressureModel.extraWellSolOutput = true;
seqModel.transportModel.extraWellSolOutput = true;


timer = tic();
[ws_split, states_split, report_split] = simulateScheduleAD(state, ...
   seqModel, schedule, 'NonLinearSolver', solver);
t_split = toc(timer);



%% Run fully implicit simulation

cprsolver = CPRSolverAD('ellipticSolver', amgSolver);

timer = tic();
[ws_fi, states_fi, report_fi] = simulateScheduleAD(state, model, ...
   schedule, 'LinearSolver', cprsolver);
t_fi = toc(timer);



%% Plot the well solutions and simulator states
% We setup interactive viewers for both well solutions and the reservoir
% states.

time = {report_fi.ReservoirTime, report_split.ReservoirTime}; 
ws = {ws_fi, ws_split};

plotWellSols(ws, time, 'datasetnames', {'FI', 'sequential'})


%% Plot each solution

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


%% Plot saturation difference

figure; plotToolbar(G, states_fi{1}.s - states_split{1}.s)


%% Plot flux field

for i = 1:numel(states_split)
    states_split{i}.v = faceFlux2cellVelocity(G, states_split{i}.flux(:, 1));
end

figure;
plotToolbar(G, states_split);
axis tight off;
colorbar




