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


%% Plot grid

figure;
plotCellData(G, rock.perm(:,1)./(milli*darcy), 'facealpha', 0.7);
view([-16,20]); axis tight; colorbar;
plotWell(G, schedule.control(1).W);
xlabel('x-axis');ylabel('y-axis');zlabel('z-axis');



%% Load polymer fluid from example in ad-blackoil

fn    = 'ad-blackoil/examples/polymerExamples/POLYMER.DATA';
deck  = readEclipseDeck(fn);
deck  = convertDeckUnits(deck);
fluid = initDeckADIFluid(deck);

%model = OilWaterPolymerModel(G, rock, fluid, 'inputdata', deck);
model = TwoPhaseOilWaterModel(G, rock, fluid, 'inputdata', deck);

% Remove gass component from wells
% schedule.control.W(1).compi = [1 0];
% schedule.control.W(2).compi = [0 1];

% Test
%schedule.control.W(2).lims.bhp = 300*barsa;
% schedule.control.W(2).type = 'bhp';
% schedule.control.W(2).val  = 300*barsa;
% 
% ijk = gridLogicalIndices(G);
% schedule.control.W(1).cells = find(ijk{1}==1 & ijk{2}==1);
% schedule.control.W(1).dir   = repmat('Z',3,1);
% schedule.control.W(1).cstatus = ones(3,1);
% schedule.control.W(2).cells = find(ijk{1}==max(ijk{1}) & ...
%     ijk{2}==max(ijk{2}));
% schedule.control.W(2).dir   = repmat('Z',3,1);
% schedule.control.W(2).cstatus = ones(3,1);

W = addWell([], G, rock, find(ijk{1}==1 & ijk{2}==1), 'type', 'rate', ...
    'val', 300*barsa, 'radius', 0.01, 'dir', 'z', 'sign', 1, ...
    'comp_i', [1 0], 'name', 'INJECTOR');
W = addWell(W, G, rock, find(ijk{1}==max(ijk{1}) & ...
    ijk{2}==max(ijk{2})), 'type', 'rate', ...
    'val', 4000*meter^3/day, 'radius', 0.01, 'dir', 'z', 'sign', 1, ...
    'comp_i', [0 0], 'name', 'INJECTOR');


% Create schdedule
dt = 1*day*ones(10,1);
schedule = simpleSchedule(dt, 'Wells', W);

% schedule.step.control = schedule.step.control(1:nsteps);
% schedule.step.val     = schedule.step.val(1:nsteps);

% Inject polymer
%schedule.control.W(1).poly = 0;

% Set initial state
state.s = state.s(:,1:2);
%state.c = zeros(G.cells.num, 1);
%state.cmax = zeros(G.cells.num, 1);

% Oil rel-perm from 2p OW system.
% Needed by equation implementation function 'eqsfiOWExplictWells'.
%fluid.krO = fluid.krOW;



%% TEST: Pure Oil/Water case

% %model = TwoPhaseOilWaterModel(G, rock, fluid, 'inputdata', deck);
% model = OilWaterPolymerModel(G, rock, fluid, 'inputdata', deck);
% 
% % Injector: water injection
% schedule.control.W(1).compi = [1 0];
% schedule.control.W(1).type  = 'rate';
% schedule.control.W(1).val   = 3000 * (meter^3) / day;
% 
% % Producer
% schedule.control.W(2).compi = [0 1];
% schedule.control.W(2).type  = 'bhp';
% schedule.control.W(2).val   = 300*barsa;
% 
% % Set initial state
% state.s = state.s(:,1:2);
% state.c = zeros(G.cells.num, 1);
% state.cmax = zeros(G.cells.num, 1);
% 
% % Reduce schedule
% nsteps = 10;
% schedule.step.control = schedule.step.control(1:nsteps);
% schedule.step.val     = schedule.step.val(1:nsteps);




%% Plot initial data

if false
    %%
    figure;
    plotCellData(G, state.s(:,1), 'facealpha', 0.7);
    %plotCellData(G, state.pressure./barsa, 'facealpha', 0.7);
    view([-16,20]); axis tight; colorbar;
    plotWell(G, schedule.control(1).W);
    xlabel('x-axis');ylabel('y-axis');zlabel('z-axis');
end

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

return;

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

for i = 2:2
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




