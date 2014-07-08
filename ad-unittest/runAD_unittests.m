mrstModule add ad-unittest

%% Run SPE 1
tests = matlab.unittest.TestSuite.fromClass(?TestSPE1);
result_spe1 = run(tests);

%% Run SPE 3
tests = matlab.unittest.TestSuite.fromClass(?TestSPE3);
result_spe3 = run(tests);

%% Run SPE 9
tests = matlab.unittest.TestSuite.fromClass(?TestSPE9);
result_spe9 = run(tests);

%% Run oil/water test
tests = matlab.unittest.TestSuite.fromClass(?TestSimpleOW);
result_ow = run(tests);

%% Run oil/water test
tests = matlab.unittest.TestSuite.fromClass(?TestEGG);
result_egg = run(tests);

%% Run tests on the nonlinear solver
tests = matlab.unittest.TestSuite.fromClass(?Test_simulateScheduleAD);
result_sched = run(tests);

%%
tests = matlab.unittest.TestSuite.fromClass(?ResultHandlerTest);
res = run(tests);

