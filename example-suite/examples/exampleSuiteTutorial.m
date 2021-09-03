mrstModule add example-suite ad-core ad-props ad-blackoil mrst-gui
mrstVerbose on

%% List examples
listExampleSuite();

%% Get and plot example
name = 'qfs_wo'; % Choose a name from the table
example = MRSTExample(name);
example.plot(example.model.rock);

%% Simulate
problem = example.getPackedSimulationProblem();
simulatePackedProblem(problem);

%% Interactive plotting
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
example.plot(states);

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
