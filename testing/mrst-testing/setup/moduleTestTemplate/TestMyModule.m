classdef TestMyModule < matlab.unittest.TestCase
%TestMyModule  Template unit test class for MRST modules.
%
% USAGE:
%   Copy this file into your module's ``tests/`` directory and rename both
%   the file and the class to match your module (e.g. ``TestMyModule.m``).
%   Add one method per test in the ``methods (Test)`` block.
%
% CONVENTIONS:
%   - File must be named  Test*.m  or  *Test.m  (e.g. TestMyModule.m).
%   - Place in  <moduleRoot>/tests/  so getUnitTestSuiteMRST discovers it.
%   - Each test method must be in the ``methods (Test)`` block.
%   - Use ``tc.assert*`` / ``tc.verify*`` for assertions.
%   - Use ``tc.assumeTrue(condition, msg)`` to skip gracefully when
%     optional dependencies are missing.
%
% SEE ALSO:
%   getUnitTestSuiteMRST, runPreReleaseTests

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

    % ---------------------------------------------------------------
    % Run once before all tests in this class
    % ---------------------------------------------------------------
    methods (TestClassSetup)
        function loadModules(tc) %#ok<MANU>
            % Add required MRST modules
            % mrstModule add ad-core my-module
        end
    end

    % ---------------------------------------------------------------
    % Run before each test method
    % ---------------------------------------------------------------
    methods (TestMethodSetup)
        function setupTest(tc) %#ok<MANU>
            % Per-test setup, e.g. create a small test grid
        end
    end

    % ---------------------------------------------------------------
    % Test methods
    % ---------------------------------------------------------------
    methods (Test)

        function examplePassingTest(tc)
            % A trivial sanity check to verify the test infrastructure works.
            tc.verifyEqual(1 + 1, 2);
        end

        function exampleAssertionTest(tc)
            % Demonstrate a numerical assertion with tolerance.
            expected = pi;
            computed = 4 * atan(1);
            tc.verifyEqual(computed, expected, 'AbsTol', 1e-14);
        end

        function exampleSkipTest(tc)
            % Demonstrate skipping when an optional feature is unavailable.
            hasPCT = ~isempty(ver('parallel'));
            tc.assumeTrue(hasPCT, ...
                'Parallel Computing Toolbox not available; skipping.');
            % ... rest of test ...
        end

        function exampleErrorTest(tc)
            % Verify that a function throws an expected error.
            tc.verifyError(@() error('myid:myerror', 'oops'), 'myid:myerror');
        end

    end
end
