mrstModule add ad-fi deckformat mrst-gui ad-core ad-blackoil

% Because several examples use the SPE1 dataset, the initial setup is
% delegated to a helper function. See the inside for documentation.
[G, rock, fluid, deck, state] = setupSPE1();

% Determine the model automatically from the deck. It should be a
% three-phase black oil model with gas dissoluton.
model = selectModelFromDeck(G, rock, fluid, deck);

model %#ok, intentional display

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(G, model, rock, deck);
%% Run base case 
[wellSols, states, report] = simulateScheduleAD(state, model, schedule);

%% Make a smaller schedule
% Because the SPE1 benchmark only has a single well configuration during
% the entire simulation, we are (relatively) free to choose timesteps. To
% demonstrate this, we create a simpler schedule consisting of a single
% very long timestep.

schedule_small = schedule;
schedule_small.step.val     = sum(schedule.step.val);
schedule_small.step.control = 1;

%% Run various number of target iterations
% Small step to get the solver started
rampup = 1*day;
targetIts = [4 8 15 25];
[reports, ws] = deal(cell(numel(targetIts), 1));
for i = 1:numel(targetIts)
    % Set up time selection class
    timestepper = ...
       IterationCountTimeStepSelector('targetIterationCount', targetIts(i),...
                                      'minRelativeAdjustment', sqrt(eps),...
                                      'maxRelativeAdjustment', inf, ...
                                      'firstRampupStep',       rampup, ...
                                      'verbose', true);

    % Instansiate a non-linear solver with the timestep class as a
    % construction argument.
    nonlinear = NonLinearSolver('timeStepSelector', timestepper, ...
                                'maxiterations', 100);
    % Solve and store results.
    [ws{i}, ~, reports{i}] = simulateScheduleAD(state, model, schedule_small,...
                        'nonlinearSolver', nonlinear, 'outputMinisteps', true);
end

%% Combine output
% We now have reports and solutions for both base case and the cases using
% automatic timesteping. Combine the results into 
l = arrayfun(@(x) ['Target: ', num2str(x), ' its'], targetIts, 'UniformOutput', false);
l = horzcat('Basecase', l);

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
    plot(time{i}, timesteps{i}, '--o', 'linewidth', 2, 'color', c(i, :))
    grid on
end
legend(l)

%% Compare the solutions interactively
plotWellSols(wsols, time, 'datasetnames', l)

%% Find the number of iterations and simulation time taken for all cases
% Since the timesteps produce substeps we have to find report->control step
% reports.
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
