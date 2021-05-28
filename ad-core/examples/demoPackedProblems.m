%% Managing simulations: Restart, packed problems and more
% When using MRST to simulate larger or many problems, it is best to ensure
% that the results are stored to disk in a predictible manner. Otherwise,
% output from a long simulation may be lost if Matlab is closed for
% whatever reason.
mrstModule add ad-core ad-blackoil ad-props mrst-gui
%% Set up a simple problem with quadratic and linear relperm
% We set up a simple 1D displacement problem for purposes of this
% demonstration. We create the same scenario with two different fluid
% models: A version with linear relative permeability and one with
% quadratic.
name = '1d';
G = cartGrid([100, 1], [1000, 1000]);
rock = makeRock(G, 0.1*darcy, 0.3);
G = computeGeometry(G);
fluid_1 = initSimpleADIFluid('mu', [1, 1]*centi*poise, 'n', [1, 1], 'phases', 'wo', 'rho', [1000, 500]);

model_1 = TwoPhaseOilWaterModel(G, rock, fluid_1);
time = 10*year;
irate = sum(model_1.operators.pv)/time;
W = [];
W = addWell(W, G, rock, 1, 'type', 'rate', 'val', irate, 'comp_i', [1, 0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 100*barsa, 'comp_i', [1, 0]);

dt = rampupTimesteps(time, time/100);
schedule = simpleSchedule(dt, 'W', W);

state0 = initResSol(G, 150*barsa, [0, 1]);

% Second fluid model
fluid_2 = initSimpleADIFluid('mu', [1, 1]*centi*poise, 'n', [2, 2], 'phases', 'wo', 'rho', [1000, 500]);
model_2 = TwoPhaseOilWaterModel(G, rock, fluid_2);
%% Define two simulations
% We create two atomic packed simulation problems. This is similar to using
% simulateScheduleAD, with the addition of a few parameters describing the
% problem.

% Define base name for all problems in this script
BaseName = 'test_packed';
% Define linear flux
problem_1 = packSimulationProblem(state0, model_1, schedule, BaseName, ...
                                'Name', 'linear_relperm', ...
                                'Description', '1D displacement with linear flux');
% Define nonlinaer flux
problem_2 = packSimulationProblem(state0, model_2, schedule, BaseName, ...
                                'Name', 'quadratic_relperm', ...
                                'Description', '1D displacement with nonlinear flux');
%% Run simulation
% This will only simulate what is not yet simulated, and can restart aborted
% simulations in a seamless manner by using the ResultHandler class
% internally.
problems = {problem_1, problem_2};
[ok, status] = simulatePackedProblem(problems);
%% We can get output, just as if we simulated in the standard manner
[ws_1, states_1, reports_1] = getPackedSimulatorOutput(problem_1);
[ws_2, states_2, reports_2] = getPackedSimulatorOutput(problem_2);
%% For many packed problems, we can also get the data in a single call
[all_ws, all_states, all_reports, names, report_time] = getMultiplePackedSimulatorOutputs(problems);
%% Plot well results
plotWellSols(all_ws, report_time, 'DatasetNames', names, 'field', 'qWs', 'SelectedWells', 2);
%% Remove data (with prompts)
clearPackedSimulatorOutput(problems)
%% We can also launch simulations in the background
% This is useful for running multiple cases. This functionality does not
% require any Matlab add-ons, since it creates an additional Matlab session
% which runs in the background. Please note that the additional sessions
% will terminate upon completion on error, but for larger cases this may
% take some time. You should not launch more sessions than your workstation
% can handle!
%
% First, remove data (without prompt before deletion)
clearPackedSimulatorOutput(problems, 'prompt', false)
% Launch simulations
for i = 1:numel(problems)
    info = simulatePackedProblemBackground(problems{i}, 'verbose', true);
end
% Monitor the simulations running in the background
monitorBackgroundSimulations(problems);

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
