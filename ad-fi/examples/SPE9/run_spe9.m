%% Ninth Comparative Solution Project
% This example runs the model from Killough, J. E. 1995. Ninth SPE
% comparative solution project: A reexamination of black-oil simulation. In
% SPE Reservoir Simulation Symposium,  12-15 February 1995, San Antonio,
% Texas. SPE 29110-MS, doi: 10.2118/29110-MS
mrstModule add ad-fi ad-core deckformat mrst-gui ad-props
mrstVerbose true

%% Read and process file.
% This <ninth SPE comparative solution project> consists of a water
% injection problem in a highly heterogenous reservoir. There is one
% injector and 25 producers. The problem is set up to be solved using a
% black-oil model. The data set we provide is a modified version of input
% files belonging to the <http://www.ntnu.edu/studies/courses/TPG4535
% course in reservoir engineering and petrophysics> at NTNU (Trondheim,
% Norway) and available at
% <http://www.ipt.ntnu.no/~kleppe/pub/SPE-COMPARATIVE/ECLIPSE_DATA/>.

%% Read input files
current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'SPE9.DATA');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

fluid = initDeckADIFluid(deck);
gravity on

%% Initial solution
p0  = deck.SOLUTION.PRESSURE;
sw0 = deck.SOLUTION.SWAT;
sg0 = deck.SOLUTION.SGAS;
s0  = [sw0, 1-sw0-sg0, sg0];
rs0 = deck.SOLUTION.RS;
rv0 = 0;

state = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);
clear k p0 s0 rs0

%% Set up the schedule
% There are three control periods and the controls change after 300 and 360
% days. However, well control can switch if the well constraints are not
% fulfilled (this is the default option. If you do not want this, set the
% option 'allowControlSwitching' to false when calling initADIsystem). In
% the case of the schedule we are considering, the wells are controlled by
% rate (oil rate for the producers and water rate for the injectors). There
% is a minimum value for bottom hole pressure for the producers (here, 68.9
% bar for all producing wells). It means that if the bottom hole pressure
% (bhp) which is needed to obtain the given rate is below the minimum bhp,
% then the well changes control type and becomes controlled by pressure. For
% the injectors, we have a maximum value for the bhp (here, 275.79 bar).

schedule = deck.SCHEDULE;
system = initADISystem(deck, G, rock, fluid, 'cpr', true);
system.nonlinear.cprBlockInvert = false;
% convergence is overall better for quite strict limits on update
system.stepOptions.drsMax = .2;
system.stepOptions.dpMax  = .2;
system.stepOptions.dsMax  = .05;
% gmres tol needs to be quite strict
system.nonlinear.cprRelTol = 1e-3;

%% Plot reservoir and wells
figure, 
plotCellData(G, log10(convertTo(rock.perm(:,1), milli*darcy))); view(3)
h = colorbar('horiz');
set(h, 'XTickLabel', num2str(10 .^ (get(h, 'XTick') .')), ...
    'YTick', 0.5, 'YTickLabel', '[mD]', ...
    'Position', [0.13, 0.05, 0.77, 0.03]);

W = processWells(G,rock, schedule.control(1));
nm = [ { 'I' }, ...
       arrayfun(@(w) sprintf('P-%02d', w - 1), 2 : numel(W), ...
                'UniformOutput', false) ];
[ W.name ] = nm{:}; clear nm

plotWell(G, W(1),    'color', 'b', 'LineWidth', 3, 'FontSize', 14);
plotWell(G, W(2:end),'color', 'k', 'LineWidth', 3, 'FontSize', 14);
axis tight, view(40,55);

%% Run schedule
% We use the fully implicit solver
timer = tic;
[wellSols, states, iter] = ...
   runScheduleADI(state, G, rock, system, schedule);
toc(timer)

%% Plot solutions
% We opt for a simple volume plot of the oil saturation.
%

figure(1)
clf
view(35, 40);
for i = 2:numel(states)
   [az, el] = view();
   clf
   plotWell(G, W, 'FontSize', 10, 'LineWidth', 1);
   plotCellData(G, states{i}.s(:,2));
   time = sum(schedule.step.val(1:i-1));
   title(['Oil saturation. Step ', num2str(i), ...
          ' (', formatTimeRange(time) ')'])
   axis tight off
   view(az, el);
   pause(0.1)
end

%% Plot of the production and injection rates.
% We observe that the rates are not constant in the control periods given
% initially ([0,300], [300,360], [360,900] days). For some timesteps, the
% well controls have switched from rate to bhp.

T     = convertTo(cumsum(schedule.step.val), day);

% Put the well solution data into a format more suitable for plotting
[qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);
injInx = 1;
prdInx = 2:26;

figure(1)
plot(T, convertTo(qWs(:, injInx), meter^3/day), '-*')
xlabel('Time [days]'), title('Water injection rate [m^3/day]')

figure(2)
plot(T, convertTo(bhp(:, prdInx), barsa), '-*'),
xlabel('Time [days]'), title('BHPs [bar]')

figure(3)
for i = 1:25
   subplot(5,5,i)
   plot(T, convertTo(qOs(:, prdInx(i)), meter^3/day), '-*');
   xlabel('Time [days]');
   title(sprintf('Oil production rate (%s) [m^3/day]', ...
                 W(prdInx(i)).name));
end
set(gcf, 'position', [100, 100, 1500, 1000])

figure(4)
plot(T, convertTo(qGs(:, prdInx), 1e3*meter^3/day), '-*')
xlabel('Time [days]'), title('Gas production rates [10^3 m^3/day]')
a = axis; axis([a(1:2), -15, a(4)])
line(a(1:2), [0, 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)

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
