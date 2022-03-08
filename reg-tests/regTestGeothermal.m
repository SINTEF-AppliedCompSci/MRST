%% Regression tests for geothermal module
% This script sets ups and runs a set of regression tests for the
% geothermal module. Useful to track changes and/or bugs introduced during
% development.

%% Add modules
mrstModule add ad-core ad-props ad-blackoil compositional
mrstModule add geothermal
mrstModule add linearsolvers
mrstModule add test-suite
mrstModule add upr
mrstVerbose off

%% Set up regression tests
% Define test cases for regression testing
tests = {
    % Test case name        % Test case options
    'qfs_geothermal'      , {'Name', 'qfs'};
    'qfs_geothermal'      , {'NaCl', true, 'Name', 'qfs_nacl'};
    'small_egs_geothermal', {};
    'htates_geothermal'   , {};
};
% Set up regression test group
rtest = RegressionTestGroup('geothermal', tests(:,1), 'regTestOpt', tests(:,2));

%% Run regression tests
rtest.runRegressionTests();

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>