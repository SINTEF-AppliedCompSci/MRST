function rts = getModuleRegressionTests()
%getModuleRegressionTests  Template: return regression tests for this module.
%
% DESCRIPTION:
%   Place this file at the root of your MRST module to register regression
%   tests that will be picked up by getAllRegressionTests and run as part
%   of Tier 2 pre-release testing.
%
%   Each element of the returned cell array must be a RegressionTest
%   object.  The RegressionTest constructor accepts either a test-case name
%   (string matching a setup function in ad-unittest/setup-functions) or a
%   TestCase object.
%
% RETURNS:
%   rts - Cell array of RegressionTest objects.
%
% EXAMPLE:
%   To add a regression test for the ''egg_wo'' test case:
%
%       rts = { RegressionTest('egg_wo', 'group', 'mymodule') };
%
% SEE ALSO:
%   RegressionTest, TestCase, getAllRegressionTests, runPreReleaseTests

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

    % Replace the example below with your module's actual test cases.
    rts = {};

    % Example (uncomment and adapt):
    % rts = {
    %     RegressionTest('my_test_case', 'group', 'mymodule'), ...
    % };
end
