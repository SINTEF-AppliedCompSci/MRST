mrstModule add geothermal
mrstModule add test-suite
mrstModule add ad-core ad-props ad-blackoil
mrstModule add compositional
mrstModule add mrst-gui

mrstVerbose on

%%
cs = 'd';
test = TestCase('benchmark_1d_geothermal', 'case', cs);
test.model.thermalFormulation = 'enthalpy';

%%
nls = NonLinearSolver();
nls.useRelaxation = true;
nls.maxIterations = 500;
nls.maxTimestepCuts = 20;
problem = test.getPackedSimulationProblem('NonLinearSolver', nls, 'Name', cs);
simulatePackedProblem(problem, 'restartStep', 1);

%%
close all
[ws, st, rep] = getPackedSimulatorOutput(problem);
time = cumsum(test.schedule.step.val);
[~, ix] = min(abs(test.options.tplot - time));
ix = min(ix, numel(st));
test.plot(st{ix}, 'field', 'T')
test.plot(st, 'field', 'T')

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