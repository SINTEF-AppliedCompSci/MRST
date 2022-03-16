%% Add required modules
mrstModule add test-suite
mrstModule add ad-core ad-props ad-blackoil compositional
mrstModule add upr
mrstModule add geothermal
mrstModule add mrst-gui
mrstVerbose on
checkHashSettings();

%% Load test case
test = TestCase('small_egs_geothermal'); test.plot();

%% Plot setup
test.figure();
plotGrid(test.model.G, test.model.G.cells.tag, 'faceColor', [1,1,1]*0.8, 'edgeColor', 'none');
plotGrid(test.model.G, 'faceColor', 'none', 'edgeAlpha', 0.1);
test.plotWells(); test.setAxisProperties(gca);
camlight(); axis off

%% Simulate problem
problem = test.getPackedSimulationProblem('useHash', true);
simulatePackedProblem(problem, 'restartStep', 1);

%% Simulate with different rock properties
test2 = test;
% Increase rock thermal conductivity by a factor 5
test2.model.rock.lambdaR = test2.model.rock.lambdaR*5;
% Update operators
test2.model = test2.model.setupOperators();
% Simulate problem
problem2 = test2.getPackedSimulationProblem('useHash', true);
simulatePackedProblem(problem2, 'restartStep', 1);

%% Plot solutions
[wellSols , states , reports ] = getPackedSimulatorOutput(problem);
[wellSols2, states2, reports2] = getPackedSimulatorOutput(problem2);
test.plot(states);
test.plot(states2);

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
