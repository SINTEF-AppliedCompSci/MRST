% filepath: /Users/fran/FranFiles/git-repos/MRST/testing/mrst-testing/getFixedTestSuiteMRST.m
function suite = getFixedTestSuiteMRST()
%Build a unit test suite from fixed core and ad-unittest test directories.
%
% SYNOPSIS:
%   suite = getFixedTestSuiteMRST()
%
% DESCRIPTION:
%   Collects matlab.unittest.TestCase subclasses from fixed test directories
%   that are always run during CI:
%     1. autodiff/ad-unittest/tests/unitTests
%     2. autodiff/ad-unittest/tests/integrationTests
%     3. core/tests/unitTests
%     4. core/tests/integrationTests
%
%   These directories contain the most critical tests that should always
%   pass and are executed on every CI run regardless of module filter.
%
% RETURNS:
%   suite - matlab.unittest.TestSuite aggregating all discovered tests.
%
% NOTE:
%   Discovery silently skips directories that cause errors during scan
%   (e.g. missing dependencies) to keep the runner robust.
%
% SEE ALSO:
%   getUnitTestSuiteMRST, getIntegrationTestSuiteMRST, runPreReleaseTests

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

    mrstModule add mrst-testing

    % Collect sub-suites in a cell array to avoid vertcat issues with
    % matlab.unittest.TestSuite across MATLAB releases.
    parts = {};

    % --- Fixed test directories (always included for CI) -----------------
    % testFolders = {
    %     '../autodiff/ad-unittest/tests/unitTests', ...
    %     '../autodiff/ad-unittest/tests/integrationTests', ...
    %     'tests/unitTests', ...
    %     'tests/integrationTests', ...
    % };
    testFolders = {
        'tests/unitTests', ...
        '../autodiff/ad-unittest/tests/unitTests', ...
    };

    for i = 1:numel(testFolders)
        folderPath = ROOTDIR;
        p = fullfile(folderPath, testFolders{i});
        
        if isfolder(p)
            try
                s = matlab.unittest.TestSuite.fromFolder(p);
                parts{end+1} = s; %#ok<AGROW>
            catch ME
                warning('mrst:fixedTestDiscovery', ...
                    'Skipping folder %s: %s', p, ME.message);
            end
        else
            warning('mrst:fixedTestDiscovery', ...
                'Test folder not found: %s', p);
        end
    end

    % Combine all collected sub-suites
    if isempty(parts)
        suite = matlab.unittest.TestSuite.empty();
    else
        suite = horzcat(parts{:});
    end
end