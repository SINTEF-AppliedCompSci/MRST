mrstModule add ad-unittest ad-core ad-blackoil

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
tests = matlab.unittest.TestSuite.fromClass(?TestSimulateScheduleAD);
result_sched = run(tests);

%%
clear
tests = matlab.unittest.TestSuite.fromClass(?ResultHandlerTest);
res = run(tests);

%%
tests = matlab.unittest.TestSuite.fromClass(?TestAdjoints);
result_adjoints = run(tests);

%% Run all simulation tests
folder = mrstPath('query', 'ad-unittest');

tests = matlab.unittest.TestSuite.fromFolder(fullfile(folder, 'test_sim'));
results = run(tests);