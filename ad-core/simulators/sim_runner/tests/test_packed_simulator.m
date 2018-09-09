mrstModule add ad-core ad-blackoil ad-props mrst-gui
%% Set up a simple problem with quadratic and linear relperm
name = '1d';
G = cartGrid([100, 1], [1000, 1000]);
rock = makeRock(G, 0.1*darcy, 0.3);
G = computeGeometry(G);
fluid_1 = initSimpleADIFluid('mu', [1, 1]*centi*poise, 'n', [1, 1], 'phases', 'wo', 'rho', [1000, 500]);

model_1 = TwoPhaseOilWaterModel(G, rock, fluid_1);
time = 10*year;
irate = sum(model_1.operators.pv)/time;
W = [];
W = addWell(W, G, rock, 1, 'type', 'rate', 'val', irate, 'comp_i', [1, 0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 100*barsa, 'comp_i', [1, 0]);

dt = rampupTimesteps(time, time/100);
schedule = simpleSchedule(dt, 'W', W);

state0 = initResSol(G, 150*barsa, [0, 1]);

% Second fluid model
fluid_2 = initSimpleADIFluid('mu', [1, 1]*centi*poise, 'n', [2, 2], 'phases', 'wo', 'rho', [1000, 500]);
model_2 = TwoPhaseOilWaterModel(G, rock, fluid_2);
%% Define two simulations
% Both linear and nonlinear relperms for the same basic configuration
problem_1 = packSimulationProblem(state0, model_1, schedule, 'test_packed', ...
                                'Name', 'linear_relperm');
                            
problem_2 = packSimulationProblem(state0, model_2, schedule, 'test_packed', ...
                                'Name', 'quadratic_relperm');
%% Run simulation
% This will only simulate what is needed, and can restart aborted
% simulations seemlessly using the ResultHandler class
problems = {problem_1, problem_2};
[ok, status] = simulatePackedProblem(problems);
%% We can get output, just as if we simulated in the standard manner
[ws_1, states_1, reports_1] = getPackedSimulatorOutput(problem_1);
[ws_2, states_2, reports_2] = getPackedSimulatorOutput(problem_2);
%% Remove data (with prompts)
clearPackedSimulatorOutput(problems)