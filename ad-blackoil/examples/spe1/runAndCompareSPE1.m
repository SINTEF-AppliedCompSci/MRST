mrstModule add ad-fi ad-props deckformat mrst-gui ad-core ad-blackoil

% Because several examples use the SPE1 dataset, the initial setup is
% delegated to a helper function. See the inside for documentation.
[G, rock, fluid, deck, state] = setupSPE1();

% Determine the model automatically from the deck. It should be a
% three-phase black oil model with gas dissoluton.
model = selectModelFromDeck(G, rock, fluid, deck);

model %#ok, intentional display

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(model, deck);

%% Plot well and permeability
% The permeability consists of three layers going from high to low
% permeability along the z axis. The wells are completed in the upper and
% lower layer for the injector and producer respectively. To get a well
% object, we simply process the first control from the schedule.
%
% Note that a schedule is not necessary to solve problems using the fully
% implicit solver: solvefiADI is capable of taking a well object directly
% and solving for a single time step in a manner similar to the other MRST
% solvers.
figure;

% Pick the only well control present
W = schedule.control(1).W;

plotCellData(G, convertTo(rock.perm(:,1), milli*darcy), ...
             'FaceAlpha', 0.5, 'EdgeAlpha', 0.3, 'EdgeColor', 'k');
plotWell(G, W);
title('Permeability (mD)')
axis tight;
view(35, 40);
colorbar('SouthOutside');

%% Run the entire schedule
[wellSols, states, report] = simulateScheduleAD(state, model, schedule);

%% Plot the well solutions and simulator states
% We setup interactive viewers for both well solutions and the reservoir
% states.

plotWellSols(wellSols, report.ReservoirTime)

figure;
plotToolbar(G, states)
plotWell(G, schedule.control(1).W)
axis tight
view(-10, 60)

%% Make plots to compare solution with a commercial simulator
% Load summary from binary file and find indices of the producer and
% injector.
load SPE1_smry

inj = find([wellSols{1}.sign] == 1);
prod = find([wellSols{1}.sign] == -1);

% Since there are zero values in the first step of the summary, we ignore
% the first entry to get better plot axes.
ind = 2:118;

% Put the well solution data into a format more suitable for plotting
[qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);

% Get timesteps for both the reference and the MRST run
T = convertTo(cumsum(schedule.step.val), year);
Tcomp =  smry.get(':+:+:+:+', 'YEARS', ind);


%% Plot Producer Gas/Oil ratio
% The most interesting part of the SPE1 case is the gas/oil ratio at the
% producer. We convert the field units and plot the dimensionless ratio.
% As should be apparent from the figure, the implicit solver is able to
% qualitatively reproduce the same outflow profile.
figure(3)
clf
ecl = convertFrom(smry.get('PRODUCER', 'WGOR', ind), 1000*ft^3/stb)';
mrst = qGs(:,prod)./qOs(:,prod);

hold on
plot(T, mrst)
plot(Tcomp, ecl, 'r');
legend({'MRST', 'Eclipse'})
xlabel('Time [Years]')
title('Gas rate / Oil rate')

%% Plot Injector Bottom Hole Pressure
% The injector is rate controlled and so the bottom hole pressure is solved
% in the implicit loop. Plot it to verify accuracy.
figure(4)
clf
ecl = convertFrom(smry.get('PRODUCER', 'WBHP', ind), psia)';
mrst = bhp(:,prod);
hold on
plot(T,     convertTo(mrst, barsa))
plot(Tcomp, convertTo(ecl, barsa), 'r');
legend({'MRST', 'Eclipse'})
xlabel('Time [Years]')
ylabel('bar')
title('Bottom hole pressure (Producer)')

%% Plot Injector Bottom Hole Pressure
figure(5)
clf
ecl = convertFrom(smry.get('INJECTOR', 'WBHP', ind), psia)';
mrst = bhp(:,inj);
hold on
plot(T,     convertTo(mrst, barsa))
plot(Tcomp, convertTo(ecl, barsa), 'r');
legend({'MRST', 'Eclipse'})
xlabel('Time [Years]')
ylabel('bar')
title('Bottom hole pressure (Injector)')
