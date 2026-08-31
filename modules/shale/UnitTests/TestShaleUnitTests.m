classdef TestShaleUnitTests < matlab.unittest.TestCase
%TestShaleUnitTests  Wraps the shale module's legacy unit test functions
%  as a matlab.unittest.TestCase so they are auto-discovered by
%  getUnitTestSuiteMRST.
%
% Each test method calls the corresponding ``test*.m`` function, which
% returns 1 for pass and 0 for fail.  The result is surfaced as a proper
% matlab.unittest assertion.
%
% NOTE:
%   The shale tests require the hfm, compositional, upr and mrst-gui
%   modules as well as optional toolboxes (Statistics, ADFNE).  Tests that
%   cannot run due to missing dependencies are skipped with assumeTrue so
%   they do not count as failures.

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
        function loadShaleModule(tc) %#ok<MANU>
            % Ensure the shale module (and its test helpers) are on the path
            mrstModule add shale
            modPath = mrstPath('query', 'shale');
            if ~isempty(modPath)
                addpath(fullfile(modPath, 'UnitTests'));
            end
        end
    end

    methods (Test)
        function testEDFMNFRWrapper(tc)
            tc.runLegacyTest(@testEDFMNFR, 'testEDFMNFR');
        end

        function testpEDFMNFRWrapper(tc)
            tc.runLegacyTest(@testpEDFMNFR, 'testpEDFMNFR');
        end

        function testPEDFM_EDFM_explicitWrapper(tc)
            tc.runLegacyTest(@testPEDFM_EDFM_explicit, 'testPEDFM_EDFM_explicit');
        end

        function testDiffusionWrapper(tc)
            tc.runLegacyTest(@testDiffusion, 'testDiffusion');
        end

        function testSorptionWrapper(tc)
            tc.runLegacyTest(@testSorption, 'testSorption');
        end

        function testGangiWrapper(tc)
            tc.runLegacyTest(@testGangi, 'testGangi');
        end
    end

    methods (Access = private)
        function runLegacyTest(tc, fn, name)
            % Verify the function is available on the path
            tc.assumeTrue(~isempty(which(func2str(fn))), ...
                sprintf('Test function ''%s'' not found. Check shale module path.', name));
            try
                pass = fn();
            catch ME
                % Surface runtime errors as test failures with a useful message
                tc.assertFail(sprintf('Test ''%s'' threw an error: %s', name, ME.message));
                return
            end
            tc.verifyEqual(pass, 1, ...
                sprintf('Legacy test ''%s'' returned failure (pass=%d).', name, pass));
        end
    end
end
