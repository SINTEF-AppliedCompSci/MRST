mrstModule add ad-core ad-blackoil deckformat ad-props linearsolvers example-suite
mrstVerbose on
if ~exist('name', 'var')
    name = 'spe1';
end
switch name
    case 'spe1'
        [G, rock, fluid, deck, state] = setupSPE1();
    case 'spe9'
        [G, rock, fluid, deck, state] = setupSPE9();
    otherwise
        error('Case %d not supported', name);
end
[state0, model, schedule, nls] = initEclipseProblemAD(deck, 'TimestepStrategy', 'none');
%% Fully-implicit
% The default discretization in MRST is fully-implicit. Consequently, we
% can use the model as-is.
implicit = packSimulationProblem(state0, model, schedule, name, ...
    'NonLinearSolver', nls, 'Name', 'Fully-Implicit');
%% Explicit solver
% This solver has a time-step restriction based on the CFL condition in
% each cell. The solver estimates the time-step before each solution.
model_explicit = setTimeDiscretization(model, 'explicit');
explicit = packSimulationProblem(state0, model_explicit, schedule, name, ...
    'NonLinearSolver', nls, 'Name', 'Explicit');
%% Adaptive implicit
% We can solve some cells implicitly and some cells explicitly based on the
% local CFL conditions. For many cases, this amounts to an explicit
% treatment far away from wells or driving forces. The values for
% estimated composition CFL and saturation CFL to trigger a switch to
% implicit status can be adjusted.
model_aim = setTimeDiscretization(model, 'adaptive-implicit');
aim = packSimulationProblem(state0, model_aim, schedule, name,...
    'NonLinearSolver', nls, 'Name', 'Adaptive-Implicit (AIM)');
%% Simulate the problems
problems = {implicit, explicit, aim};
simulatePackedProblem(problems);
%% Get output and plot the well results
% There are oscillations in the well curves. Increasing the CFL limit
% beyond unity will eventually lead to oscillations.
[ws, states, reports, names, T] = getMultiplePackedSimulatorOutputs(problems);
plotWellSols(ws, T, 'datasetnames', names);
%% Plot the fraction of cells which are implicit
imp_frac = cellfun(@(x) sum(x.implicit)/numel(x.implicit), states{3});
figure;
plot(100*imp_frac);
title('Number of implicit cells in AIM')
ylabel('Percentage of implicit cells')
xlabel('Step index')
%% Plot the time-steps taken
% The implicit solvers use the control steps while the explicit solver must
% solve many substeps.
figure; hold on
markers = {'o', '.', 'x'};
for i = 1:numel(reports)
    dt = getReportMinisteps(reports{i});
    x = (1:numel(dt))/numel(dt);
    plot(x, dt/day, 'marker', markers{i})
end
xlabel('Timestep length [day]')
legend(names)

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
