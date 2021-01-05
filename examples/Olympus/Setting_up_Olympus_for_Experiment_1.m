 
mrstModule add example-suite ad-core ad-props ad-blackoil mrst-gui
mrstVerbose on

%% Get and plot example
name = 'olympus_field_wo'; % Choose a name from the table
example = MRSTExample(name);
%example.plot(example.model.rock);

%% Simulate
problem = example.getPackedSimulationProblem();
simulatePackedProblem(problem);

%% Interactive plotting
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
