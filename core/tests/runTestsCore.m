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
%
% DEPRECATED:
%   This function is deprecated. Use runPreReleaseTests.m instead, which 
%   provides the new three-tier testing framework:
%     - Tier 1: Fast unit tests (< 1 min each)
%     - Tier 2: Integration / regression tests (< 10 min each)
%     - Tier 3: Example smoke tests
%
%   NEW USAGE:
%     % Run all core tests (tiers 1 and 2)
%     runPreReleaseTests('tier', [1 2], 'modules', {'ad-core'})
%
%     % Run only tier 1
%     runPreReleaseTests('tier', 1, 'modules', {'ad-core'})
%
%     % Run in parallel
%     runPreReleaseTests('tier', 1, 'modules', {'ad-core'}, 'parallel', true)

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

    % Parse options
    opt = struct('tier', 'all', ...
                 'runInParallel', false, ...
                 'stopOnFailure', false);
    opt = merge_options(opt, varargin{:});

    % Show deprecation warning
    warning('mrst:deprecatedFunction', ...
        ['runTestsCore is deprecated. Use runPreReleaseTests instead.\n', ...
         '  See runTestsCore.m for migration examples.']);

    % Convert legacy tier names to new tier numbers
    tiers = [1 2];  % Default to all tiers (was 'all')
    if strcmpi(opt.tier, 'tier1')
        tiers = 1;
    elseif strcmpi(opt.tier, 'tier2')
        tiers = 2;
    end

    % Call the new runPreReleaseTests function
    results = runPreReleaseTests(...
        'tier', tiers, ...
        'modules', {'ad-core'}, ...
        'parallel', opt.runInParallel, ...
        'stopOnFail', opt.stopOnFailure);
end

