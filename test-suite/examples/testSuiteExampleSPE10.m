%% Add required modules
mrstModule add test-suite
mrstModule add spe10
mrstModule add ad-core ad-props ad-blackoil
mrstModule add mrst-gui
mrstVerbose on
checkHashSettings()

%% Load the entire model
test = TestCase('spe10_wo'); test.plot();

%% Load subsets of the model
test = TestCase('spe10_wo', 'layers', 'tarbert'   ); test.plot(); % Tarbert
test = TestCase('spe10_wo', 'layers', 'upper_ness'); test.plot(); % Upper Ness
test = TestCase('spe10_wo', 'layers', 13          ); test.plot(); % Layer 13

%% Run the test
% We run the test for layer 13 using a nonlinear solver with line search
nls = NonLinearSolver('useLineSearch', true);
problem = test.getPackedSimulationProblem('NonLinearSolver', nls, 'useHash', true);
simulatePackedProblem(problem, 'restartStep', 1);

%% Run a modified version of the test
% We then modify the test by shutting in one of the producers and lowering
% the injection rate by 50 %. The TestCase class will automatically detect
% that this is adifferent setup than the one above, and store the
% simulation results in a different folder.
test2 = test;                                  % Copy the test
test2.schedule.control(1).W(2).status = false; % Shut in well P2
test2.schedule.control(1).W(5).val ...         % Reduce injection rate
    = test2.schedule.control(1).W(5).val*0.5;
problem2 = test2.getPackedSimulationProblem('NonLinearSolver', nls, 'useHash', true);
simulatePackedProblem(problem2, 'restartStep', 1);

%% Load and compare results
[~, states ] = getPackedSimulatorOutput(problem );
test.plot(states); colormap(bone);
[~, states2] = getPackedSimulatorOutput(problem2);
test2.plot(states2); colormap(bone);

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
