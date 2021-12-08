%% Add required modules
mrstModule add example-suite
mrstModule add spe10
mrstModule add ad-core ad-props ad-blackoil
mrstModule add mrst-gui
mrstVerbose on

%% Load the entire model
test = TestCase('spe10_wo'); test.plot();

%% Load subsets of the model
test = TestCase('spe10_wo', 'layers', 'tarbert'   ); test.plot(); % Tarbert
test = TestCase('spe10_wo', 'layers', 'upper_ness'); test.plot(); % Upper Ness
test = TestCase('spe10_wo', 'layers', 13          ); test.plot(); % Layer 13

%% Run the test
% We run the test for layer 13 using a nonlinear solver with line search
nls = NonLinearSolver('useLineSearch', true);
problem = test.getPackedSimulationProblem('NonLinearSolver', nls);
simulatePackedProblem(problem);

%% Run a modified version of the test
% We then modify the test by shutting in one of the producers and lowering
% the injection rate by 50 %. The TestCase class will automatically detect
% that this is adifferent setup than the one above, and store the
% simulation results in a different folder.
test2 = test;                                  % Copy the test
test2.schedule.control(1).W(2).status = false; % Shut in well P2
test2.schedule.control(1).W(5).val ...         % Reduce injection rate
    = test2.schedule.control(1).W(5).val*0.5;
problem2 = test2.getPackedSimulationProblem('NonLinearSolver', nls);
simulatePackedProblem(problem2);

%% Load and compare results
[~, states ] = getPackedSimulatorOutput(problem );
test.plot(states); colormap(bone);
[~, states2] = getPackedSimulatorOutput(problem2);
test2.plot(states2); colormap(bone);