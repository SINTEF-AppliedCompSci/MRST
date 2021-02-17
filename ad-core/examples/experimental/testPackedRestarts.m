mrstModule add ad-core ad-blackoil ad-props mrst-gui
%% Set up a simple problem with quadratic and linear relperm
G = cartGrid([100, 1], [1000, 1000]);
rock = makeRock(G, 0.1*darcy, 0.3);
G = computeGeometry(G);
fluid = initSimpleADIFluid('mu', [1, 1]*centi*poise, 'n', [1, 1], 'phases', 'wo', 'rho', [1000, 500]);

model = TwoPhaseOilWaterModel(G, rock, fluid);
time = 10*year;
irate = sum(model.operators.pv)/time;
W = [];
W = addWell(W, G, rock, 1, 'type', 'rate', 'val', irate, 'comp_i', [1, 0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 100*barsa, 'comp_i', [1, 0]);

dt = time/100;
timesteps = rampupTimesteps(time/2, dt, 0);
schedule = simpleSchedule(timesteps, 'W', W);

state0 = initResSol(G, 150*barsa, [0, 1]);
%% Define two simulations
problem = packSimulationProblem(state0, model, schedule, 'test_packed_restart');


nls_mini = NonLinearSolver('timestepselector', SimpleTimeStepSelector('maxTimestep', dt/4));
problem_restart = packSimulationProblem(state0, model, schedule, 'test_packed_restart',...
                                    'NonLinearSolver', nls_mini, ...
                                    'ExtraArguments', {'outputMinisteps', true});
%% First, test regular problems
clearPackedSimulatorOutput(problem, 'prompt', false);
%% Test running from the start
simulatePackedProblem(problem)
%% Check when problem is done
simulatePackedProblem(problem)
%% Test restarts after missing results
clearPackedSimulatorOutput(problem, 'start', 30, 'prompt', false);
simulatePackedProblem(problem)
%% Test specific restart
simulatePackedProblem(problem, 'restartStep', 10)


%% Next, repeat test for ministep output
clearPackedSimulatorOutput(problem_restart, 'prompt', false);
%% Test running from the start
simulatePackedProblem(problem_restart)
%% Check when problem is done
simulatePackedProblem(problem_restart)
%% Test restarts after missing results
clearPackedSimulatorOutput(problem_restart, 'start', 30*4 + 2, 'prompt', false);
simulatePackedProblem(problem_restart)
%% Test specific restart
simulatePackedProblem(problem_restart, 'restartStep', 10)

%%
[~, states] = getPackedSimulatorOutput(problem_restart);
t = cellfun(@(x) x.time, states);
figure(1); clf
plot(t)