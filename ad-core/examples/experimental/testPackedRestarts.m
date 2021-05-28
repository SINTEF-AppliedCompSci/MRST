mrstModule add ad-core ad-blackoil ad-props mrst-gui
%% Set up a simple problem with quadratic and linear relperm
G = cartGrid([100, 1], [1000, 1000]);
rock = makeRock(G, 0.1*darcy, 0.3);
G = computeGeometry(G);
fluid = initSimpleADIFluid('mu', [1, 1]*centi*poise, 'n', [1, 1], 'phases', 'wo', 'rho', [1000, 500]);

model = TwoPhaseOilWaterModel(G, rock, fluid);
time = 10*year;
irate = sum(model.operators.pv)/time;
W = [];
W = addWell(W, G, rock, 1, 'type', 'rate', 'val', irate, 'comp_i', [1, 0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 100*barsa, 'comp_i', [1, 0]);

dt = time/100;
timesteps = rampupTimesteps(time/2, dt, 0);
schedule = simpleSchedule(timesteps, 'W', W);

state0 = initResSol(G, 150*barsa, [0, 1]);
%% Define two simulations
problem = packSimulationProblem(state0, model, schedule, 'test_packed_restart');


nls_mini = NonLinearSolver('timestepselector', SimpleTimeStepSelector('maxTimestep', dt/4));
problem_restart = packSimulationProblem(state0, model, schedule, 'test_packed_restart',...
                                    'NonLinearSolver', nls_mini, ...
                                    'ExtraArguments', {'outputMinisteps', true});
%% First, test regular problems
clearPackedSimulatorOutput(problem, 'prompt', false);
%% Test running from the start
simulatePackedProblem(problem)
%% Check when problem is done
simulatePackedProblem(problem)
%% Test restarts after missing results
clearPackedSimulatorOutput(problem, 'start', 30, 'prompt', false);
simulatePackedProblem(problem)
%% Test specific restart
simulatePackedProblem(problem, 'restartStep', 10)


%% Next, repeat test for ministep output
clearPackedSimulatorOutput(problem_restart, 'prompt', false);
%% Test running from the start
simulatePackedProblem(problem_restart)
%% Check when problem is done
simulatePackedProblem(problem_restart)
%% Test restarts after missing results
clearPackedSimulatorOutput(problem_restart, 'start', 30*4 + 2, 'prompt', false);
simulatePackedProblem(problem_restart)
%% Test specific restart
simulatePackedProblem(problem_restart, 'restartStep', 10)

%%
[~, states] = getPackedSimulatorOutput(problem_restart);
t = cellfun(@(x) x.time, states);
figure(1); clf
plot(t)
%%
[~, states] = getPackedSimulatorOutput(problem_restart);
%%
problems = {problem, problem_restart};
[ws, states, reports, names, T, grids, act] = getMultiplePackedSimulatorOutputs(problems, 'minSteps', 20);

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
