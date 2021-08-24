%% Automated Time-Step Control
% The ad-core module offers several options for automatic time-step
% control. The simplest approach is to run the time steps as prescribed and
% use a heuristic algorithm to reduce the time step if the number of
% nonlinear iterations exceeds a prescribed limit; this is commonly
% referred to as the Appleyard chop. In the modified version, one also
% monitors the saturation increments and reduces the time step if the
% saturation changes more than a prescribed limit within a single time
% step.
%
% MRST also offers more a sophisticated time-step control that monitors the
% convergence of the nonlinear solves and adjusts the time-step so that the
% number of nonlinear iterations stays close to a prescribed target. In
% this example, we will use the SPE1 benchmark to demonstrate this time
% step control.

mrstModule add ad-props deckformat mrst-gui ad-core ad-blackoil example-suite

%% Load and run the SPE1 model
% Because several examples use the SPE1 dataset, the initial setup is
% delegated to a helper function. See the inside for documentation.
[G, rock, fluid, deck, state] = setupSPE1();

% Determine the model automatically from the deck. It should be a
% three-phase black oil model with gas dissoluton.
model = selectModelFromDeck(G, rock, fluid, deck);

model %#ok, intentional display

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(model, deck);

% Run base case 
[wellSols, states, report] = simulateScheduleAD(state, model, schedule);

%% Time-step control with iteration targets
% Because the SPE1 benchmark uses the same well configuration throughout
% the entire simulation, we are (relatively) free to choose timesteps. To
% demonstrate this, we create a simpler schedule consisting of a single
% very long timestep and then let the time-step control automatically
% select the length of the individual time steps so that the number of
% iterations for each step stays close to a prescribed target. We will run
% four simulations with target values 4, 8, 15, and 25.
%
% Use the compressSchedule routine to reduce the schedule to the minimal
% number of timesteps required to honor the distinct control steps. In this
% case, we only have a single control step. To get a stable solution, we
% also need to start with a small time step to get the solver started.
% Here, we start with a time step of one day and then gradually increase
% the length of the time step to stay close to our iteration target.

% Reduce to a minimal schedule
schedule_small = compressSchedule(schedule);

% Run all four simulations
rampup = 1*day;
targetIts = [4 8 15 25];
[reports, ws] = deal(cell(numel(targetIts), 1));
for i = 1:numel(targetIts)
    % Set up time selection class
    timestepper = ...
       IterationCountTimeStepSelector('targetIterationCount', targetIts(i),...
                                      'minRelativeAdjustment', sqrt(eps),...
                                      'maxRelativeAdjustment', 4, ...
                                      'firstRampupStep',       rampup, ...
                                      'verbose', true);

    % Instantiate a nonlinear solver with the timestep class as a
    % construction argument.
    nonlinear = NonLinearSolver('timeStepSelector', timestepper, ...
                                'maxiterations', 4*targetIts(i));
    % Solve and store results.
    [ws{i}, ~, reports{i}] = simulateScheduleAD(state, model, schedule_small,...
                        'nonlinearSolver', nonlinear, 'outputMinisteps', true);
end

%% Combine output
% We now have reports and solutions for the base case and the cases using
% automatic time-stepping. We extract the length of each time step from the
% simulation reports and plot versus time
l = arrayfun(@(x) ['Target: ', num2str(x), ' its'], targetIts, 'UniformOutput', false);
l = horzcat('Base case', l);

wsols = vertcat({wellSols}, ws);

timesteps = cell(numel(reports) + 1, 1);
timesteps{1} = schedule.step.val;
for i = 1:numel(reports)
    % Read out the actually used timesteps
    [~, t] = convertReportToSchedule(reports{i}, schedule_small);
    timesteps{i+1} = t;
end
% Sum up all timesteps to get time at each datapoint
time = cellfun(@cumsum, timesteps, 'UniformOutput', false);

% Plot the timestep lengths
figure;
hold on
c = lines(numel(timesteps));
for i = 1:numel(timesteps);
    plot(time{i}/day, timesteps{i}/day, '-o', 'linewidth', 2, 'color', c(i, :))
    grid on
end
legend(l,'Location','NorthWest')

%% Compare solution accuracy (well curves)
% Looking at the various well curves, we see that whereas the water rate
% and water production is accurately resolved by all simulations, the
% larger time steps fail to accurately resolve the strong gradients in the
% oil and gas production as the gas front approaches the producer. On the
% other hand, using a target of 15 or 25 iterations allows time steps
% between one and two years, which are extremely, and one should generally
% not expect that this will give very accurate solutions.
plotWellSols(wsols, time, 'datasetnames', l);

%% Compare computational efficiency
% To get an overall view of the computational cost of each simulation, we
% need to extract the number of iterations and computational time of each
% substep. These can be extracted from the report->control step reports.
reps = vertcat({report}, reports);

[iterations, timesteps] = deal(zeros(numel(reps), 1));
for i = 1:numel(reps)
    r = reps{i};
    timesteps(i) = sum(r.SimulationTime);
    for j = 1:numel(reps{i})
        for k = 1:numel(r.ControlstepReports)
            if r.Converged
                iterations(i) = iterations(i) + r.ControlstepReports{k}.Iterations;
            end
        end
    end
end
figure;
bar([iterations, timesteps])
legend('Non-linear iterations', 'Time taken')
set(gca, 'XTicklabel', l)

%% Conclusion
% The automatic time-step control enables us to run very large time steps.
% However, since this will also introduce large errors in the temporal
% discretization, a reasonable compromise between accuracy and efficiency
% would most likely be to select an iteration target somewhere in the
% interval 3 to 8.

%% Copyright notice

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
