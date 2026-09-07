function suite = getCoreUnitTestSuiteMRST()
%Build the matlab.unittest.TestSuite for MRST's base (core/) unit tests
%
% SYNOPSIS:
%   suite = getCoreUnitTestSuiteMRST()
%
% DESCRIPTION:
%   Returns the suite of `matlab.unittest`-conformant tests covering
%   `core/` functionality. Unlike the analogous autodiff suite builders
%   (see `getUnitTestSuiteMRST`), this does not blindly scan `core/` with
%   `TestSuite.fromFolder`, since `core/` also contains files that merely
%   look like tests by name (e.g. old demo/plotting scripts) but are not
%   `matlab.unittest`-conformant, which would break folder-based discovery.
%   Instead, known-good test files are listed explicitly below. Add new
%   `core/` unit test files to this list as they are written -- see
%   TESTING.md for how to add one.
%
% RETURNS:
%   suite - A matlab.unittest.TestSuite covering all listed core test files.
%
% SEE ALSO:
%   `getUnitTestSuiteMRST`, `getIntegrationTestSuiteMRST`, `runMRSTTests`

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

    root = ROOTDIR();  % NOTE: ROOTDIR() returns the 'core/' directory itself
    testFiles = {
        fullfile(root, 'utils', 'gridtools', 'tests', 'gridToolsTest.m'), ...
        fullfile(root, 'utils', 'equil', 'simpleEquilibriumTest.m'), ...
    };

    import matlab.unittest.TestSuite;
    suite = [];
    for i = 1:numel(testFiles)
        suite = [suite, TestSuite.fromFile(testFiles{i})]; %#ok<AGROW>
    end
end
