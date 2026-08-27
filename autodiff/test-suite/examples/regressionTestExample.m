%% Add modules
mrstModule add test-suite
mrstModule add mrst-testing
mrstModule add ad-core ad-props ad-blackoil
mrstVerbose off % Supress extensive command window output

%% Set up a single regression test
% Regression tests are set up using MRSTRegressionTest.  The setup function
% can be given directly and any setup arguments are passed via setupArgs.
rtest = MRSTRegressionTest(@gravity_segregation, ...
                           'group', 'examples', ...
                           'name',  'gravity_segregation', ...
                           'setupArgs', {'ncells', 21*3});
disp(rtest);
% You can also point the test at a packed-problem builder directly.
rtest2 = MRSTRegressionTest(@gravity_segregation, ...
                            'group', 'examples', ...
                            'name',  'gravity_segregation_custom', ...
                            'setupArgs', {'ncells', 21*3});

%% Running a single regression test
% Running the test is done with a one-line command. This will first check
% for existing results for the same test. If no results exist, the function
% will issue a warning saying that the test will be inconclusive. After the
% test has been run and compared to existing results (if any), existing
% results from before the test are deleted. If you want to keep results
% from before the test, you can set the property `deletePrevious` to false.
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
% You can define a regression test using specific packed-problem inputs,
% such as a non-standard nonlinear solver, with the optional input argument
% `problemInput`. These are passed directly to `packSimulationProblem`.
nls  = NonLinearSolver('useLineSearch', true);
rtest = MRSTRegressionTest(@qfs_wo, ...
                           'group', 'examples', ...
                           'name',  'qfs_wo', ...
                           'setupArgs', {'dt', 100*day}, ...
                           'problemInput', {'NonLinearSolver', nls});
report = rtest.runRegressionTest(); %#ok

%% Defining a group of regression tests
% In larger projects, it is convenient to keep a cell array of
% MRSTRegressionTest objects in getModuleRegressionTests.m.
rtest1 = MRSTRegressionTest(@qfs_wo, 'group', 'qfs-tests', ...
                            'name', 'qfs_wo_1', 'setupArgs', {'nkr', 1});
rtest2 = MRSTRegressionTest(@qfs_wo, 'group', 'qfs-tests', ...
                            'name', 'qfs_wo_2', 'setupArgs', {'nkr', 2});
rtest3 = MRSTRegressionTest(@qfs_wo, 'group', 'qfs-tests', ...
                            'name', 'qfs_wo_3', 'setupArgs', {'nkr', 3});
report = {rtest1.runRegressionTest(), rtest2.runRegressionTest(), rtest3.runRegressionTest()};

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.
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
