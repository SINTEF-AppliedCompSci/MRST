function skipList = getSkippedExamples()
%getSkippedExamples  Template: return example basenames to skip for this module.
%
% DESCRIPTION:
%   Place this file at the root of your MRST module to declare which
%   examples should be excluded from automated smoke testing (Tier 3).
%
%   The returned cell array should contain example file basenames (without
%   path, with or without the .m extension).  These examples are filtered
%   out by MRSTExampleTests when building the test suite.
%
%   Prefer per-file MRST_TEST_OPTIONS annotations for fine-grained control:
%
%       %{
%       MRST_TEST_OPTIONS
%       interactive: true    % skip in automated runs
%       tags:        slow    % excluded by default
%       %}
%
%   Use this function only for batch exclusions that apply to many files or
%   cannot be annotated (e.g. third-party files you do not own).
%
% RETURNS:
%   skipList - Cell array of example basenames (strings).
%
% SEE ALSO:
%   mrstExampleOptions, MRSTExampleTests, runPreReleaseTests

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

    skipList = {
        % Add basenames of examples to skip, e.g.:
        % 'myVeryLongExample',   % reason: takes > 30 min
        % 'myInteractiveDemo',   % reason: requires user input
    };
end
