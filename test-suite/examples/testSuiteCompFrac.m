%% Minimal example showing how to run multiple examples using MRSTExample

%% Add required modules
mrstModule add test-suite
mrstModule add ad-core ad-props compositional
mrstModule add mrst-gui
mrstVerbose on

%%
test = TestCase('fractures_compositional'); test.plot(); colormap(pink);

%%
problem = test.getPackedSimulationProblem();
simulatePackedProblem(problem);

%%
test2 = test;
test2.schedule.control(1).W(1).val = test2.schedule.control(1).W(1).val*0.25;
problem2 = test2.getPackedSimulationProblem();
simulatePackedProblem(problem2);

%%
[wellSols, states, reports] = getPackedSimulatorOutput(problem); test.plot(states);
[wellSols2, states2, reports2] = getPackedSimulatorOutput(problem2); test2.plot(states2);