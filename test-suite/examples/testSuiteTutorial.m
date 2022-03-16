%% Test suite tutorial
% In this tutorial, we will go through the basic concepts of a test suite,
% including how to call a test case setup function, how to identify a text
% case by it's hash value, and how to run and visualize a test case using
% the TestCase class.

%% Add necessary modules
mrstModule add test-suite
mrstModule add ad-core ad-props ad-blackoil
mrstModule add mrst-gui
mrstVerbose off
checkHashSettings();

%% Using a test case setup function
% We consider a quarter five-spot setup posed on a 1000 x 1000 m domain.
% This is implemented in the test case setup function `qfs_wo`. Providing
% the optional input argument `fullSetup` determines if we should return
% the options and descriptions only (false), or the full setup (true). See
% the function `testcase_template` to get an overwiew of the general body
% of a test case setup function.
setup0 = qfs_wo(false, 'nkr', 1); % Get options and depscription only
disp(setup0);

setup  = qfs_wo(true, 'nkr', 1);  %#ok Get full setup
setup  = qfs_wo('nkr', 1);        % Omitting fullSetup also gives full setup
disp(setup);

%% Setting up a test case
% The TestCase class implements useful functionality for setting up,
% running and visualizing a test case. Each test case can be uniquely
% determined by a hash value, which makes it easy to distinguish the
% instances of the same test case set up with different options. We get the
% test case hash using the method `getTestCaseHash`. This takes an optional
% argument `fullSetup` that determines if we should include the entire test
% case in the hash computation, or just the name, description and options.
% QFS with linear relative permeabilities
test = TestCase('qfs_wo', 'nkr', 1);
hash = test.getTestCaseHash(false); disp(hash)
hash = test.getTestCaseHash();      disp(hash)
% QFS with quadratic relative permeabilities
test = TestCase('qfs_wo', 'nkr', 2);
hash = test.getTestCaseHash(false); disp(hash)
hash = test.getTestCaseHash();      disp(hash)
% Change permeability in first cell
test.model.rock.perm(1) = test.model.rock.perm(1)*2;
hash = test.getTestCaseHash(false); disp(hash)
hash = test.getTestCaseHash();      disp(hash)

%% Simulating a test case
% This block of code contains all the necessary commands for setting up and
% running a test case.
test = TestCase('qfs_wo', 'nkr', 2, 'ncells', 75); % Quadr relperms, [75,75] cells
problem = test.getPackedSimulationProblem();       % Get test case problem
simulatePackedProblem(problem, 'restartStep', 1);  % Simulate

%% Visualizing the test case
% The TestCase method plot wraps around plotToolbar and generates visually
% pleasing plots of the data.
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
test.plot(states); colormap(bone);
% You can also use the TestCase methods to generate visually pleasing plots
% without using plotToolbar
test.figure();
plotCellData(test.model.G, states{1}.pressure, 'edgeAlpha', 0.2);
test.plotWells();
test.setAxisProperties(gca);

%% Saving and loading a test case to/from disk
% The TestCase class offers functionality for saving test cases to disk,
% and loading them at a later time. This is useful when working with test
% cases that takes a long time to set up. By default, a test case is stored
% to disk in a subfolder of mrstDataDirectory() using the test case hash
% computed from the name, description and options as a file name. We can
% also specify a different directory and name through optional input
% arguments.
test.save('directory', []   , ... % Empty gives the default directory
          'name'     , []   , ... % Empty gives the test case hash
          'prompt'   , true);     % Setting this to true (default) will 
                                  % bring up a prompt with directory, file
                                  % name and file size before saving
                                  
% Next time we attempt to set up the same test case, the class will first
% look for a saved version and load this. With verbose = true, the user is
% informed whether the test case was found on disk, and how long it took to
% load.
test = TestCase('qfs_wo', 'nkr', 2, 'ncells', 75, 'verbose', true);

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
