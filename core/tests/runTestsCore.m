function results = runTestsCore(varargin)
%Run tiered tests for MRST core functionality.
%
% SYNOPSIS:
%   results = runTestsCore()
%   results = runTestsCore('tier', 'tier1')
%
% OPTIONAL PARAMETERS:
%   tier          - 'all' (default), 'tier1', or 'tier2'.
%   runInParallel - Whether to execute the suite in parallel. Default false.
%   stopOnFailure - Whether to stop after the first failure. Default false.

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

   opt = struct('tier', 'all', ...
                'runInParallel', false, ...
                'stopOnFailure', false);
   opt = merge_options(opt, varargin{:});

   import matlab.unittest.TestRunner;
   import matlab.unittest.TestSuite;
   import matlab.unittest.plugins.StopOnFailuresPlugin;
   import matlab.unittest.selectors.HasTag;

   testRoot = fileparts(mfilename('fullpath'));
   suite    = TestSuite.fromFolder(testRoot);

   if ~strcmpi(opt.tier, 'all')
      suite = suite.selectIf(HasTag(lower(opt.tier)));
   end

   runner = TestRunner.withTextOutput;
   if opt.stopOnFailure
      runner.addPlugin(StopOnFailuresPlugin);
   end

   if opt.runInParallel
      results = runner.runInParallel(suite);
   else
      results = runner.run(suite);
   end
end

