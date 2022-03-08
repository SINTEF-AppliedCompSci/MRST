%% Add required modules
mrstModule add test-suite
mrstModule add ad-core ad-props ad-blackoil compositional
mrstModule add upr
mrstModule add geothermal
mrstModule add mrst-gui
mrstVerbose on
checkHashSettings();

%% Load test case
test = TestCase('small_egs_geothermal'); test.plot();

%% Plot setup
test.figure();
plotGrid(test.model.G, test.model.G.cells.tag, 'faceColor', [1,1,1]*0.8, 'edgeColor', 'none');
plotGrid(test.model.G, 'faceColor', 'none', 'edgeAlpha', 0.1);
test.plotWells(); test.setAxisProperties(gca);
camlight(); axis off

%% Simulate problem
problem = test.getPackedSimulationProblem('useHash', true);
simulatePackedProblem(problem);

%% Simulate with different rock properties
test2 = test;
% Increase rock thermal conductivity by a factor 5
test2.model.rock.lambdaR = test2.model.rock.lambdaR*5;
% Update operators
test2.model = test2.model.setupOperators();
% Simulate problem
problem2 = test2.getPackedSimulationProblem('useHash', true);
simulatePackedProblem(problem2);

%% Plot solutions
[wellSols , states , reports ] = getPackedSimulatorOutput(problem);
[wellSols2, states2, reports2] = getPackedSimulatorOutput(problem2);
test.plot(states);
test.plot(states2);