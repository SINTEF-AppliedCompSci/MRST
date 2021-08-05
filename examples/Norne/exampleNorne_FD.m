mrstModule add example-suite ad-core ad-props ad-blackoil mrst-gui
mrstModule add fd-publications kpn
mrstVerbose on


%% Get and plot example
name = 'norne_simple_wo';%'norne_fd_publication'; % Choose a name from the table
example = MRSTExample(name);


%% Simulate
problem = example.getPackedSimulationProblem();
simulatePackedProblem(problem);

%% Interactive plotting
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
example.plot(states);

