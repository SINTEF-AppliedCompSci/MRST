%% Script for running tests using github actions
% This file must be run in the directory battmoDir()/Tests

% Display matlab version
disp(version)

%% Display git commit ids and other stats
testdir = pwd;
[~, res] = system('git rev-parse --short HEAD');
fprintf('%s %s', pwd, res);
cd(mrstPath('ad-core'));
[~, res] = system('git rev-parse --short HEAD');
fprintf('%s', res);
cd(testdir)

global MRST_BATCH
MRST_BATCH = true;

import matlab.unittest.TestSuite;
import matlab.unittest.selectors.*;
import matlab.unittest.parameters.Parameter;
import matlab.unittest.TestRunner


%% Run tests

% Setup
mrstVerbose 'off';
stopOnError        = false;
runTestsInParallel = true;
doAssertSuccess    = true;

% Define which test cases to run
testCases = {'TestExamples'};

% Setup test suite
for itestcase = 1 : numel(testCases)

    testCase = testCases{itestcase};
    suite =  TestSuite.fromClass(meta.class.fromName(testCase));

    % if strcmp(testCase, 'TestExamples')

        % % Tests that are not supported on github
        % filenames = {'somtest'};
        % selector = HasParameter('Property', 'filename', 'Value', filenames{1});
        % for ifile = 2 : numel(filenames)
        %     filename = filenames{ifile};
        %     selector = HasParameter('Property', 'filename', 'Value', filename) | selector;
        % end
        % suite = suite.selectIf(~selector);

    % end

    suites{itestcase} = suite;

end

suite = horzcat(suites{:});

runner = testrunner('textoutput');

if stopOnError
    import matlab.unittest.plugins.StopOnFailuresPlugin;
    runner.addPlugin(StopOnFailuresPlugin);
end

% Run tests

if runTestsInParallel
    results = runner.runInParallel(suite);
else
    results = runner.run(suite);
end

t = table(results);
disp(t);

if doAssertSuccess
    assertSuccess(results);
end

