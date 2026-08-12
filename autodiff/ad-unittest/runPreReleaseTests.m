function results = runPreReleaseTests(varargin)
%Run MRST pre-release tests across one or more testing tiers.
%
% SYNOPSIS:
%   results = runPreReleaseTests()
%   results = runPreReleaseTests('tier', tiers)
%   results = runPreReleaseTests('tier', tiers, 'modules', mods)
%   results = runPreReleaseTests(..., 'excludeTags', tags)
%   results = runPreReleaseTests(..., 'writeXML', true)
%
% DESCRIPTION:
%   Three-tier pre-release testing system for MRST:
%
%     Tier 1 — Fast unit tests (< 1 min each, run on every push/PR).
%              Uses matlab.unittest.TestCase subclasses discovered via
%              getUnitTestSuiteMRST.
%
%     Tier 2 — Integration / regression tests (< 10 min each, nightly).
%              Full simulations compared against stored reference results.
%              Requires RegressionTest objects returned by modules.
%
%     Tier 3 — Example smoke tests (run before release).
%              Runs every non-interactive, non-excluded example via
%              MRSTExampleTests / getExampleIntegrationTestSuiteMRST.
%
% OPTIONAL PARAMETERS:
%   'tier'         - Which tiers to run.  Scalar or vector, e.g. 1, [1 2],
%                    or [1 2 3].  Default: [1 2 3].
%   'modules'      - Cell array of module names.  When given, only tests
%                    belonging to these modules are executed.  Default: all.
%   'excludeTags'  - Cell array of example tags to exclude in Tier 3.
%                    Default: {'slow', 'data-required'}.
%   'writeXML'     - Write JUnit XML results to ad-unittest/output/XML/.
%                    Default: false.
%   'stopOnFail'   - Stop the test run on first failure. Default: false.
%   'parallel'     - Run tests in parallel (requires Parallel Computing
%                    Toolbox). Default: false.
%
% RETURNS:
%   results - Struct array with fields:
%               tier    (double)
%               name    (char)
%               result  (matlab.unittest.TestResult array or struct)
%
% EXAMPLES:
%   % Fast pre-commit check
%   runPreReleaseTests('tier', 1);
%
%   % Nightly check
%   runPreReleaseTests('tier', [1 2]);
%
%   % Full pre-release for a specific module
%   runPreReleaseTests('tier', [1 2 3], 'modules', {'ad-blackoil'});
%
%   % Run including slow examples
%   runPreReleaseTests('tier', 3, 'excludeTags', {});
%
% SEE ALSO:
%   getUnitTestSuiteMRST, getIntegrationTestSuiteMRST,
%   getExampleIntegrationTestSuiteMRST, getAllRegressionTests,
%   mrstModuleTestStatus

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    % --- Parse options ---------------------------------------------------
    opt = struct( ...
        'tier'       , [1 2 3]              , ...
        'modules'    , {{}}                 , ...
        'excludeTags', {{'slow','data-required'}}, ...
        'writeXML'   , false                , ...
        'stopOnFail' , false                , ...
        'parallel'   , false                );

    % Allow environment variable overrides for CI usage
    opt = applyEnvOverrides(opt);
    opt = merge_options(opt, varargin{:});

    mrstModule add ad-unittest

    import matlab.unittest.TestRunner;
    import matlab.unittest.TestSuite;

    xmlDir = '';
    if opt.writeXML
        xmlDir = fullfile(mrstPath('query', 'ad-unittest'), 'output', 'XML');
        if ~isfolder(xmlDir), mkdir(xmlDir); end
    end

    results = struct('tier', {}, 'name', {}, 'result', {});

    % =====================================================================
    % TIER 1 — Unit tests
    % =====================================================================
    if any(opt.tier == 1)
        fprintf('\n%s\n Running Tier 1: Unit Tests\n%s\n', repmat('=',1,60), repmat('=',1,60));
        suite = getUnitTestSuiteMRST();
        if ~isempty(opt.modules)
            suite = filter_by_modules(suite, opt.modules);
        end
        res = run_suite(suite, 'unit', opt, xmlDir);
        results(end+1) = struct('tier', 1, 'name', 'unit', 'result', {res});
    end

    % =====================================================================
    % TIER 2 — Integration / regression tests
    % =====================================================================
    if any(opt.tier == 2)
        fprintf('\n%s\n Running Tier 2: Integration Tests\n%s\n', repmat('=',1,60), repmat('=',1,60));
        % Simulation integration tests (matlab.unittest based)
        intSuite = getIntegrationTestSuiteMRST();
        if ~isempty(opt.modules)
            intSuite = filter_by_modules(intSuite, opt.modules);
        end
        res = run_suite(intSuite, 'integration', opt, xmlDir);
        results(end+1) = struct('tier', 2, 'name', 'integration', 'result', {res});

        % RegressionTest objects
        rts = getAllRegressionTests();
        if ~isempty(rts)
            runRegressionTests(rts, opt);
        end
    end

    % =====================================================================
    % TIER 3 — Example smoke tests
    % =====================================================================
    if any(opt.tier == 3)
        fprintf('\n%s\n Running Tier 3: Example Smoke Tests\n%s\n', repmat('=',1,60), repmat('=',1,60));
        mods = opt.modules;
        if isempty(mods)
            mods = mrstPath();
        end
        [suites, suiteNames] = getExampleIntegrationTestSuiteMRST( ...
            mods, ...
            'seperateModules', true, ...
            'excludeTags'    , opt.excludeTags);
        for i = 1:numel(suites)
            if isempty(suites{i}), continue; end
            res = run_suite(suites{i}, suiteNames{i}, opt, xmlDir);
            results(end+1) = struct('tier', 3, 'name', suiteNames{i}, 'result', {res}); %#ok<AGROW>
        end
    end

    % --- Print summary ---------------------------------------------------
    printSummary(results);
end

% =========================================================================
% Helper functions
% =========================================================================

function res = run_suite(suite, name, opt, xmlDir)
    import matlab.unittest.TestRunner;
    if isempty(suite)
        fprintf('  (no tests found for ''%s'')\n', name);
        res = matlab.unittest.TestResult.empty();
        return
    end
    runner = TestRunner.withTextOutput();
    if opt.stopOnFail
        import matlab.unittest.plugins.StopOnFailuresPlugin;
        runner.addPlugin(StopOnFailuresPlugin());
    end
    if ~isempty(xmlDir)
        import matlab.unittest.plugins.XMLPlugin;
        jFile = fullfile(xmlDir, [name, '.xml']);
        runner.addPlugin(XMLPlugin.producingJUnitFormat(jFile));
    end
    if opt.parallel
        try
            res = runner.runInParallel(suite);
        catch
            warning('mrst:parallelFailed', 'Parallel run failed; retrying serially.');
            res = runner.run(suite);
        end
    else
        res = runner.run(suite);
    end
end

function runRegressionTests(rts, opt)
    fprintf('  Running %d regression test(s)...\n', numel(rts));
    for i = 1:numel(rts)
        rt = rts{i};
        try
            report = rt.runRegressionTest();
            if report.passed
                fprintf('  [PASS] %s\n', rt.name);
            else
                fprintf('  [FAIL] %s\n', rt.name);
            end
        catch ME
            fprintf('  [ERROR] %s: %s\n', rt.name, ME.message);
        end
    end
end

function suite = filter_by_modules(suite, modules)
    % Filter a TestSuite to only include tests from the given modules.
    % Works on suites that have a 'module' TestParameter.
    try
        active = arrayfun( ...
            @(x) any(strcmpi(x.Parameterization(2).Value, modules)), suite);
        suite = suite(active);
    catch
        % Suite does not have a module parameter — return as-is
    end
end

function opt = applyEnvOverrides(opt)
    % Allow CI environment variables to override defaults
    v = getenv('MRST_TEST_TIERS');
    if ~isempty(v)
        opt.tier = str2num(v); %#ok<ST2NM>
    end
    v = getenv('MRST_EXCLUDE_TAGS');
    if ~isempty(v)
        opt.excludeTags = strsplit(v, ',');
    end
    if ~isempty(getenv('MRST_WRITE_XML'))
        opt.writeXML = true;
    end
    if ~isempty(getenv('MRST_STOP_ON_FAIL'))
        opt.stopOnFail = true;
    end
end

function printSummary(results)
    fprintf('\n%s\n Pre-Release Test Summary\n%s\n', repmat('=',1,60), repmat('=',1,60));
    fprintf('%-5s  %-25s  %6s  %6s  %6s\n', 'Tier', 'Suite', 'Passed', 'Failed', 'Total');
    fprintf('%s\n', repmat('-',1,60));
    anyFailed = false;
    for i = 1:numel(results)
        r   = results(i);
        res = r.result;
        if isempty(res)
            nPass = 0; nFail = 0; nTotal = 0;
        else
            nPass  = sum([res.Passed]);
            nFail  = sum([res.Failed]);
            nTotal = numel(res);
        end
        if nFail > 0, anyFailed = true; end
        fprintf('%-5d  %-25s  %6d  %6d  %6d\n', ...
            r.tier, r.name, nPass, nFail, nTotal);
    end
    fprintf('%s\n', repmat('=',1,60));
    if anyFailed
        fprintf('RESULT: SOME TESTS FAILED\n');
    else
        fprintf('RESULT: ALL TESTS PASSED\n');
    end
end
