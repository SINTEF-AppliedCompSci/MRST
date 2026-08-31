function rts = getModuleRegressionTests()
%getModuleRegressionTests  Template: return regression tests for this module.
%
% DESCRIPTION:
%   Place this file at the root of your MRST module to register regression
%   tests that will be picked up by getAllRegressionTests and run as part
%   of Tier 2 pre-release testing.
%
%   Each element of the returned cell array must be an MRSTRegressionTest
%   object. The constructor accepts a setup function handle or a packed
%   problem builder and optional name/value arguments.
%
% RETURNS:
%   rts - Cell array of MRSTRegressionTest objects.
%
% EXAMPLE:
%   To add a regression test for the ``egg_wo`` case:
%
%       rts = { MRSTRegressionTest(@egg_wo, 'group', 'mymodule') };
%
% SEE ALSO:
%   MRSTRegressionTest, getAllRegressionTests, runPreReleaseTests

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

    rts = {};

    % Example:
    % rts = {
    %     MRSTRegressionTest(@my_regression_case, 'group', 'mymodule'), ...
    % };
end
