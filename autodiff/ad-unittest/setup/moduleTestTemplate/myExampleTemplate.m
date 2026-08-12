%% myExampleTemplate  Template for annotated MRST examples
%
% This file shows how to add a MRST_TEST_OPTIONS annotation block so that
% the pre-release testing workflow can automatically classify this example.
%
% Place the annotation block inside a ``%{ ... %}`` comment near the very
% top of the example file (within the first 50 lines).
%
%{
MRST_TEST_OPTIONS
timeout:     0            % seconds until the test is aborted; 0 = no limit
interactive: false        % set true to always skip in automated runs
tags:                     % comma-separated: slow, data-required, gpu, ...
%}
%
% Valid tag values (by convention):
%   slow           - example takes more than a few minutes; excluded by default in CI
%   data-required  - example needs external datasets not bundled with MRST
%   interactive    - example requires GUI interaction (also use interactive: true)
%   gpu            - example requires a GPU
%
% Examples:
%   A fast smoke-test example needs no annotation at all.
%
%   A long-running example:
%{
MRST_TEST_OPTIONS
tags: slow
%}
%
%   An interactive GUI demo:
%{
MRST_TEST_OPTIONS
interactive: true
%}
%
%   An example needing external data with a 5-minute timeout:
%{
MRST_TEST_OPTIONS
timeout:     300
tags:        data-required
%}

%% Setup
% (Your example code starts here)
mrstModule add ad-core

G = cartGrid([10 10], [1 1]);
G = computeGeometry(G);

disp('Template example ran successfully.');

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
