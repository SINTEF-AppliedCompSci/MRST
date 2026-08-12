function status = mrstModuleTestStatus(varargin)
%Audit pre-release test coverage for all registered MRST modules.
%
% SYNOPSIS:
%   mrstModuleTestStatus()
%   status = mrstModuleTestStatus()
%   status = mrstModuleTestStatus('verbose', true)
%
% DESCRIPTION:
%   For every registered MRST module, checks which testing artefacts are
%   present and reports a summary table to the command window.
%
%   Checked artefacts:
%     unitTests       - files matching Test*.m or *Test.m in tests/ or
%                       UnitTests/ subdirectory
%     regressionFn    - getModuleRegressionTests.m at module root
%     skipListFn      - getSkippedExamples.m at module root
%     numExamples     - total number of examples in the module
%     numAnnotated    - examples that carry a MRST_TEST_OPTIONS block
%
% OPTIONAL PARAMETERS:
%   'verbose' - Print the full table to the command window. Default: true.
%   'modules' - Cell array of module names to restrict the audit.
%               Default: all registered modules.
%
% RETURNS:
%   status - Table with one row per module and the columns described above.
%            Returned only if a return variable is requested.
%
% EXAMPLE:
%   mrstModuleTestStatus()
%   t = mrstModuleTestStatus('verbose', false);
%
% SEE ALSO:
%   runPreReleaseTests, getUnitTestSuiteMRST, getAllRegressionTests

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

    opt = struct('verbose', true, 'modules', {{}});
    opt = merge_options(opt, varargin{:});

    mods = mrstPath();
    if ~isempty(opt.modules)
        mods = intersect(mods, opt.modules, 'stable');
    end

    n = numel(mods);
    modName       = mods(:);
    unitTests     = zeros(n, 1);
    regressionFn  = false(n, 1);
    skipListFn    = false(n, 1);
    numExamples   = zeros(n, 1);
    numAnnotated  = zeros(n, 1);

    testDirs = {'tests', 'UnitTests'};

    for i = 1:n
        modPath = mrstPath('query', mods{i});
        if isempty(modPath), continue; end

        % Count unit test files
        for d = 1:numel(testDirs)
            td = fullfile(modPath, testDirs{d});
            if isfolder(td)
                f1 = dir(fullfile(td, 'Test*.m'));
                f2 = dir(fullfile(td, '*Test.m'));
                unitTests(i) = unitTests(i) + numel(f1) + numel(f2);
            end
        end

        % Regression test function
        regressionFn(i) = exist(fullfile(modPath, 'getModuleRegressionTests.m'), 'file') > 0;

        % Skip-list function
        skipListFn(i)   = exist(fullfile(modPath, 'getSkippedExamples.m'), 'file') > 0;

        % Examples
        try
            ex = mrstExamples(mods{i});
            examples = ex{1};
            numExamples(i) = numel(examples);
            ann = 0;
            for j = 1:numel(examples)
                o = mrstExampleOptions(examples{j});
                % An annotation block is present when any field differs from default
                hasAnnotation = o.interactive || o.timeout ~= 0 || ~isempty(o.tags);
                if hasAnnotation, ann = ann + 1; end
            end
            numAnnotated(i) = ann;
        catch
            % ignore errors from broken modules
        end
    end

    % Build table
    status = table(modName, unitTests, regressionFn, skipListFn, ...
                   numExamples, numAnnotated, ...
        'VariableNames', {'Module','UnitTestFiles','HasRegressionFn', ...
                          'HasSkipListFn','NumExamples','NumAnnotated'});

    if opt.verbose
        printStatusTable(status);
    end

    if nargout == 0
        clear status
    end
end

% -------------------------------------------------------------------------

function printStatusTable(t)
    fprintf('\n%s\n MRST Module Test Coverage Audit\n%s\n', ...
        repmat('=',1,72), repmat('=',1,72));
    fprintf('%-30s  %5s  %8s  %8s  %8s  %8s\n', ...
        'Module', 'UTest', 'RegFn', 'SkipFn', 'NExmpl', 'NAnnot');
    fprintf('%s\n', repmat('-',1,72));
    for i = 1:height(t)
        fprintf('%-30s  %5d  %8s  %8s  %8d  %8d\n', ...
            t.Module{i}, ...
            t.UnitTestFiles(i), ...
            tf2str(t.HasRegressionFn(i)), ...
            tf2str(t.HasSkipListFn(i)), ...
            t.NumExamples(i), ...
            t.NumAnnotated(i));
    end
    fprintf('%s\n', repmat('=',1,72));
    fprintf('Totals:  %d unit-test files  |  %d regression functions  |  %d/%d examples annotated\n', ...
        sum(t.UnitTestFiles), sum(t.HasRegressionFn), ...
        sum(t.NumAnnotated), sum(t.NumExamples));
    fprintf('%s\n\n', repmat('=',1,72));
end

function s = tf2str(v)
    if v, s = 'YES'; else, s = 'no'; end
end
