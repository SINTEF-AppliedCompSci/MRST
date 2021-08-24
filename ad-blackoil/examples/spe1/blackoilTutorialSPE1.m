%% The Odeh Benchmark (SPE1) 
% The first SPE project comparing black-oil reservoir simulators was
% organized by Odeh (1981) and describes a depletion problem with gas
% injection in a small 10x10x3 reservoir with a producer and an injector
% placed in diagonally opposite corners. The porosity is uniform and equal
% 0.3, whereas the permeability is isotropic with values 500, 50, and 200
% md in the three layers with thickness 20, 30, and 50 ft. The reservoir is
% initially undersaturated with a pressure field that is constant in each
% layer, a uniform mixture of water (Sw = 0.12) and oil (So = 0.88) with no
% initial free gas (Sg = 0.0) and a constant dissolved gas-oil ratio (Rs )
% throughout the model.
%
% The problem is specified using an industry-standard input format, and the
% solution computed by MRST is compared to that of a commercial reservoir 
% simulator (Eclipse 100). The original problem was posed to study ten
% years of production; herein, we only report and compare solutions for the
% first 1216 days.

%
% Odeh, A.S. 1981. Comparison of Solutions to a Three-Dimensional Black-Oil
% Reservoir Simulation Problem. J Pet Technol 33 (1): 13–25. SPE-9723-PA.
% http://dx.doi.org/10.2118/9723-PA 
%
mrstModule add ad-props deckformat mrst-gui ad-core ad-blackoil example-suite

%% Set up the problem
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
figure;
plotCellData(G, convertTo(rock.perm(:,1), milli*darcy), ...
             'FaceAlpha', 0.5, 'EdgeAlpha', 0.3, 'EdgeColor', 'k');
plotWell(G, schedule.control(1).W, 'radius',.5); % Pick the only well control present
title('Permeability (mD)')
axis tight, view(35, 40), colorbar('SouthOutside');

%% Run the entire schedule
% Here, we will run the schedule as it is described in the input file. Note
% that a schedule is not necessary to solve problems using the fully 
% implicit solver. The function 'solvefiADI' from the 'ad-fi' module (which
% implements fully implicit ad solvers *without* object orientation) is
% capable of taking a well structure directly and solve for a single time
% step in a manner similar e.g., to the incompressible MRST solvers.
% schedule.step.val(1) = 10*year;
nls = NonLinearSolver('useLinesearch', true);
[wellSols, states, report] = simulateScheduleAD(state, model, schedule, 'nonlinearsolver', nls);

%% Plot the well solutions and simulator states
% We setup interactive viewers for both well solutions and the reservoir
% states.
plotWellSols(wellSols, report.ReservoirTime)

figure;
plotToolbar(G, states)
plotWell(G, schedule.control(1).W)
axis tight, view(-10, 60)

%% Make plots to compare solution with a commercial simulator
% Load summary from binary file and find indices of the producer and
% injector.
load SPE1_smry

inj = find([wellSols{1}.sign] == 1);
prod = find([wellSols{1}.sign] == -1);

% Since there are zero values in the first step of the summary, we ignore
% the first entry to get better plot axes.
ind = 2:size(smry.data,2);

% Put the well solution data into a format more suitable for plotting
[qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);

% Get timesteps for both the reference and the MRST solution
T = convertTo(cumsum(schedule.step.val), year);
Tcomp =  smry.get(':+:+:+:+', 'YEARS', ind);


%% Plot producer gas/oil ratio
% The most interesting part of the SPE1 case is the gas/oil ratio at the
% producer. We convert the field units and plot the dimensionless ratio.
% As should be apparent from the figure, the implicit solver is able to
% qualitatively reproduce the same outflow profile.
figure(3)
clf
ecl = convertFrom(smry.get('PRODUCER', 'WGOR', ind), 1000*ft^3/stb)';
mrst = qGs(:,prod)./qOs(:,prod);

mrstarg = {'LineWidth', 2};
eclarg = {'ro','MarkerSize',5,'MarkerFaceColor',[.6 .6 .6]};
hold on
plot(T, mrst, mrstarg{:})
plot(Tcomp, ecl, eclarg{:});
legend({'MRST', 'Eclipse'})
xlabel('Time [Years]')
title('Gas rate / Oil rate')

%% Plot producer bottom-hole pressure
% The injector is rate controlled and so the bottom-hole pressure is solved
% in the implicit loop. Plot it to verify accuracy.
figure(4)
clf
ecl = convertFrom(smry.get('PRODUCER', 'WBHP', ind), psia)';
mrst = bhp(:,prod);
hold on
plot(T,     convertTo(mrst, barsa), mrstarg{:})
plot(Tcomp, convertTo(ecl, barsa), eclarg{:});
legend({'MRST', 'Eclipse'})
xlabel('Time [Years]')
ylabel('bar')
title('Bottom-hole pressure (producer)')

%% Plot injector bottom-hole pressure
figure(5)
clf
ecl = convertFrom(smry.get('INJECTOR', 'WBHP', ind), psia)';
mrst = bhp(:,inj);
hold on
plot(T,     convertTo(mrst, barsa), mrstarg{:})
plot(Tcomp, convertTo(ecl, barsa),  eclarg{:});
legend({'MRST', 'Eclipse'})
xlabel('Time [Years]')
ylabel('bar')
title('Bottom-hole pressure (injector)')

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
