mrstModule add example-suite ad-core ad-props ad-blackoil mrst-gui
mrstVerbose on

%% List examples
listExamples();

%% Get and plot example
name = 'qfs_wo'; % Choose a name from the table
example = MRSTExample(name);
example.plot(example.model.rock);

%% Simulate
problem = example.getPackedSimulationProblem();
simulatePackedProblem(problem);

%% Interactive plotting
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
example.plot(states);