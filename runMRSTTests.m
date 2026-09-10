function results = runMRSTTests(tier, varargin)
%Single entry point for MRST's tiered testing system
%
% SYNOPSIS:
%   results = runMRSTTests(tier)
%   results = runMRSTTests(tier, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Runs one of MRST's three test tiers and prints a command-line summary
%   of the outcome. See TESTING.md for a full overview of the tiers, how
%   they map to CI, and how to add new tests/examples.
%
% REQUIRED PARAMETERS:
%   tier - One of:
%            'unit'       - Fast, always-on tests: core/ + autodiff unit
%                            tests (autodiff/ad-unittest's test_models and
%                            test_utils suites).
%            'regression' - Heavier scenario/simulation tests that check
%                            known results stay consistent (autodiff's
%                            test_sim suite).
%            'examples'   - Runs every example script in one or more
%                            modules end-to-end as a smoke test. Requires
%                            the 'module' option.
%
% OPTIONAL PARAMETERS:
%   module         - Required for tier='examples'. A module name, a cell
%                    array of module names, or the literal string 'all' to
%                    run every registered module's examples. Not used for
%                    'unit'/'regression'.
%
%   writeXML       - If true, write JUnit-format XML test reports (one file
%                    per named sub-suite) under
%                    fullfile(mrstOutputDirectory(), 'TestOutput', ...
%                              'runMRSTTests', tier, 'XML').
%                    Default: false.
%
%   writeTAP       - If true, write TAP-format test reports similarly.
%                    Default: false.
%
%   errorOnFailure - If true, raise an error (and thus a non-zero exit code
%                    when run non-interactively) when any test failed.
%                    Default: false, so interactive use simply prints the
%                    summary without throwing. CI workflows pass true.
%
% RETURNS:
%   results - Array of matlab.unittest.TestResult for every test that ran,
%             across all sub-suites of the requested tier.
%
% EXAMPLES:
%   runMRSTTests('unit');
%   runMRSTTests('regression', 'writeXML', true);
%   runMRSTTests('examples', 'module', 'ad-blackoil');
%   runMRSTTests('examples', 'module', 'all', 'errorOnFailure', true);
%
% SEE ALSO:
%   `getCoreUnitTestSuiteMRST`, `getUnitTestSuiteMRST`,
%   `getIntegrationTestSuiteMRST`, `getExampleIntegrationTestSuiteMRST`

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

    if nargin < 1
        error('runMRSTTests:MissingTier', ...
            'Usage: runMRSTTests(tier, ...) with tier one of ''unit'', ''regression'', ''examples''.');
    end

    opt = struct(...
        'module',         {{}}, ...
        'writeXML',       false, ...
        'writeTAP',       false, ...
        'errorOnFailure', false ...
        );
    opt = merge_options(opt, varargin{:});

    % The ad-unittest suite-builder functions used below all live inside
    % the ad-unittest module's own directory, which (like any MRST module
    % under modules/, solvers/, etc.) is not on the MATLAB path until
    % explicitly activated -- do that once here so callers of this single
    % entry point never need to know about that internal detail.
    mrstModule add ad-unittest

    switch lower(tier)
        case 'unit'
            suites = { getCoreUnitTestSuiteMRST(), 'core-unit'; ...
                       getUnitTestSuiteMRST(),      'autodiff-unit' };

        case 'regression'
            suites = { getIntegrationTestSuiteMRST(), 'autodiff-regression' };

        case 'examples'
            if isempty(opt.module)
                error('runMRSTTests:MissingModule', ...
                    ['tier=''examples'' requires an explicit ''module'' ', ...
                     'option -- a module name, a cell array of module ', ...
                     'names, or ''all''. This tier is not run implicitly ', ...
                     'across every module.']);
            end
            mods = normalize_module_arg(opt.module);
            [s, n] = getExampleIntegrationTestSuiteMRST(mods, 'seperateModules', true);
            suites = [s(:), n(:)];

        otherwise
            error('runMRSTTests:UnknownTier', ...
                'Unknown tier ''%s''. Expected ''unit'', ''regression'' or ''examples''.', tier);
    end

    results = run_named_suites(suites, tier, opt);

    print_test_summary(results, suites(:, 1), suites(:, 2));

    if opt.errorOnFailure && any([results.Failed])
        error('runMRSTTests:TestsFailed', '%d of %d tests failed.', ...
              nnz([results.Failed]), numel(results));
    end
end

%--------------------------------------------------------------------------

function mods = normalize_module_arg(module)
    % Note: merge_options wraps a scalar (non-cell) override value in a
    % cell array to match the class of the 'module' option's cell-array
    % default, so a caller-supplied 'all' arrives here as {'all'}, not the
    % bare char 'all'. Handle both forms defensively.
    if ischar(module)
        module = {module};
    end
    if numel(module) == 1 && strcmpi(module{1}, 'all')
        mods = mrstPath();
    else
        mods = module;
    end
end

%--------------------------------------------------------------------------

function results = run_named_suites(suites, tier, opt)
    import matlab.unittest.TestRunner;

    if opt.writeXML
        import matlab.unittest.plugins.XMLPlugin;
        outXML = fullfile(mrstOutputDirectory(), 'TestOutput', ...
                           'runMRSTTests', lower(tier), 'XML');
        if exist(outXML, 'dir') == 0
            mkdir(outXML);
        end
    end

    if opt.writeTAP
        outTAP = fullfile(mrstOutputDirectory(), 'TestOutput', ...
                           'runMRSTTests', lower(tier), 'TAP');
        if exist(outTAP, 'dir') == 0
            mkdir(outTAP);
        end
        mrstModule add ad-unittest
    end

    results = [];
    for i = 1:size(suites, 1)
        name = suites{i, 2};
        fprintf('Running test suite %d of %d: %s\n', i, size(suites, 1), name);

        runner = TestRunner.withTextOutput;
        if opt.writeXML
            jFile = fullfile(outXML, [name, '.xml']);
            runner.addPlugin(XMLPlugin.producingJUnitFormat(jFile));
        end

        res = runner.run(suites{i, 1});
        results = [results, res]; %#ok<AGROW>

        if opt.writeTAP
            tapFile = fullfile(outTAP, [name, '.tap']);
            writeTestsTAP_YAMLISH(res, tapFile);
        end
    end
end

%--------------------------------------------------------------------------

function print_test_summary(results, suites, names)
    fprintf('\n%s\n', repmat('=', 1, 70));
    fprintf('MRST test summary\n');
    fprintf('%s\n', repmat('=', 1, 70));

    failedNames = {};
    ix = 0;
    for i = 1:numel(suites)
        n   = numel(suites{i});
        res = results(ix+1:ix+n);
        ix  = ix + n;

        nPassed  = nnz([res.Passed]);
        nFailed  = nnz([res.Failed]);
        nIncompl = nnz([res.Incomplete]);
        duration = sum([res.Duration]);

        fprintf('  %-24s %4d passed, %4d failed, %4d incomplete  (%.2fs)\n', ...
                names{i}, nPassed, nFailed, nIncompl, duration);

        fail = res([res.Failed]);
        for j = 1:numel(fail)
            failedNames{end+1} = sprintf('%s :: %s', names{i}, fail(j).Name); %#ok<AGROW>
        end
    end

    if ~isempty(failedNames)
        fprintf('\nFailed tests:\n');
        for i = 1:numel(failedNames)
            fprintf('  FAIL  %s\n', failedNames{i});
        end
    end

    nTotal    = numel(results);
    nPassed   = nnz([results.Passed]);
    nFailed   = nnz([results.Failed]);
    nIncompl  = nnz([results.Incomplete]);
    duration  = sum([results.Duration]);

    fprintf('%s\n', repmat('-', 1, 70));
    fprintf('  TOTAL: %d passed, %d failed, %d incomplete out of %d (%.2fs)\n', ...
            nPassed, nFailed, nIncompl, nTotal, duration);
    fprintf('%s\n\n', repmat('=', 1, 70));
end
