%% Minimal example showing how to run multiple examples using MRSTExample

%% Add required modules
mrstModule add example-suite
mrstModule add ad-core ad-props ad-blackoil deckformat
mrstModule add mrst-gui
mrstVerbose on

%% Simulate three SPE benchmarks
examples = {'spe1_bo', 'spe3_bo', 'spe9_bo'};
for ex = examples
    example = MRSTExample(ex{1});
    example.plot(); drawnow, pause(1);
    problem = example.getPackedSimulationProblem();
    simulatePackedProblem(problem);
end