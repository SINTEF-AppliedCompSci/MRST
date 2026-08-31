classdef TestMpfa < matlab.unittest.TestCase
%TestMpfa  Wrapper converting the fvbiot mpfaTest function-based tests to
%          matlab.unittest.TestCase so they are auto-discovered by
%          getUnitTestSuiteMRST.
%
% The test logic lives in mpfaTest (functiontests).  Here we re-run each
% named sub-test through the function-based framework and propagate any
% assertion failure as a TestCase failure.

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

    methods (TestClassSetup)
        function addFvbiotModule(tc) %#ok<MANU>
            mrstModule add fvbiot
        end
    end

    methods (Test)
        function dirichletBoundaryTest(tc)
            tc.runFunctionTest('dirichletBoundaryTest');
        end

        function mixedBoundaryTest(tc)
            tc.runFunctionTest('mixedBoundaryTest');
        end
    end

    methods (Access = private)
        function runFunctionTest(tc, testName)
            suite  = matlab.unittest.TestSuite.fromFunction(@mpfaTest);
            target = suite(arrayfun(@(s) strcmp(s.Name, testName), suite));
            tc.assertNotEmpty(target, ...
                sprintf('Sub-test ''%s'' not found in mpfaTest.', testName));
            result = run(target);
            tc.verifyTrue(result.Passed, ...
                sprintf('mpfaTest sub-test ''%s'' failed.', testName));
        end
    end
end
