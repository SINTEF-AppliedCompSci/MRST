%% Add modules
mrstModule add test-suite
mrstModule add ad-core ad-props ad-blackoil
mrstVerbose off % Supress extensive command window output
checkHashSettings();

%% Set up a single regression test
% Regression tests are set up in the same way as a TestCase, using the
% setup function name as input. Optional arguments for the test case are
% passed just as in the TestCase class.
rtest = RegressionTest('gravity_segregation', 'ncells', 21*3);
disp(rtest);
% You can also set up a regression test from an instance of the TestCase
% class. This allows for altering test case parameters before definig the
% regression test
test   = TestCase('gravity_segregation');
test.model.rock.poro(1) = 0.2;
rtest2 = RegressionTest(test); %#ok

%% Running a single regression test
% Running the test is done with a one-line command. This will first check
% for existing results for the same test. If no results exist, the function
% will issue a warning saying that the test will be inconclusive. After the
% test has been run and compared to existing results (if any), existing
% results from before the test are deleted. If you want to keep results
% from before the test, you can set the property `deleteExisting` to false.
% Results from the test are stored to disk.
report = rtest.runRegressionTest();

%% Inspect the report
% The regression test report will always have a field `passed` equal to -1
% (inconclusive, i.e., there were no existing results to compare against) 0
% (failed) or 1 (passed). If the test failed or passed, the report has a
% field `states` and `wellSols`. Each of these has a field `dvalue` where
% each field value corresponds to the accumulated deviation (in the inf
% norm) of that respective field. The `states` and `wellSols` also has a
% field `dsteps`, corresponding to the difference in number of timesteps of
% the current run to the number of steps in the previous run
disp(report)
if report.passed == 0 || report.passed == 1
    disp(report.states);
    disp(report.states.dvalues);
end

%% Running with different options
% You can define a regression test using specific inputs to
% `simulateScheduleAD`, such as a non-standard nonlinear solver, with the
% optional input argument `problemInput`. These will be passed directly to
% the TestCase method `getPackSimulationProblem`.
nls  = NonLinearSolver('useLineSearch', true);
rtest = RegressionTest('qfs_wo', 'dt', 100*day, ...
                      'problemInput', {'NonLinearSolver', nls});
report = rtest.runRegressionTest(); %#ok

%% Defining a group of regression tests
% In larger projects, it is convenient to define a set of regression tests
% that you use regularly to monitor the development process. This can be
% done using the class RegressionTestGroup. Its constructor takes in a
% group name, and a cell array of test cases (either as instances of
% TestCase, or as strings corresponding to names of test case setup
% functions), or as instances of regression tests. When an input test is
% given by a test case setup function name or an instance of TestCase,
% optional inputs can be provided with the argument `testOpt`. This should
% be a cell array with one element per test. We set up a group of
% regression tests using the quarter five-spot case with different
% Brooks-Corey relative permeability exponents.
rtest1 = TestCase('qfs_wo', 'nkr', 1, 'name', 'qfs_wo_1');
rtest2 = TestCase('qfs_wo', 'nkr', 2, 'name', 'qfs_wo_2');
rtest3 = TestCase('qfs_wo', 'nkr', 3, 'name', 'qfs_wo_3');
testGroup = RegressionTestGroup('qfs-tests', {rtest1, rtest2, rtest3});
report = testGroup.runRegressionTests();

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
