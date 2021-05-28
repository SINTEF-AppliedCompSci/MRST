%% SPE3 case using fully implicit black oil solver
% This <http://dx.doi.org/10.2118/12278-PA third SPE comparative solution
% project> consists of a gas injection problem in a small (9x9x4)
% reservoir. The problem is originally set up to be solved using a
% compositional solver. Using PVTi, the physical datas have been processed
% to obtain an equivalent blackoil problem and the resulting blackoil
% parameters are provided in the file "SPE3.DATA". The data set we provide
% is a modified version of input files belonging to the
% <http://www.ntnu.edu/studies/courses/TPG4535 course in reservoir
% engineering and petrophysics> at NTNU (Trondheim, Norway) and available
% at <http://www.ipt.ntnu.no/~kleppe/pub/SPE-COMPARATIVE/ECLIPSE_DATA/>.
% The oil can vaporize but the gas cannot dissolve in oil so that the
% gas/oil ratio remains equal to zero during the whole simulation.

%% Read input files
% The input files follow Eclipse format. MRST contains a dedicated module
% which can handle standard Eclipse keywords.

mrstModule add ad-fi ad-core deckformat ad-props

% Read and process file.
current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'SPE3.DATA');

deck = readEclipseDeck(fn);

% The deck is given in field units, MRST uses consistently the metric
% system.
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

% Create a special ADI fluid from deck, which can produce differentiated
% fluid properties.
fluid = initDeckADIFluid(deck);

% The case includes gravity
gravity on

%% Set up initial state
% The initial state corresponds to an equilibrium state between
% gravitational and capillary forces. It has been computed before and we
% load it from deck.

p0  = deck.SOLUTION.PRESSURE;
sw0 = deck.SOLUTION.SWAT;
sg0 = deck.SOLUTION.SGAS;
s0  = [sw0, 1-sw0-sg0, sg0];
rv0 = deck.SOLUTION.RV;
rs0 = 0;
state = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);   
clear k p0 s0 rv0 rs0

%% Plot wells and permeability
% The permeability is constant in each layer. There is one injecting and one
% producing well.

clf;
W = processWells(G, rock, deck.SCHEDULE.control(1));
plotCellData(G, convertTo(rock.perm(:,1), milli*darcy), ...
             'FaceAlpha', .5, 'EdgeAlpha', .3, 'EdgeColor', 'k');
plotWell(G, W, 'fontsize', 10, 'linewidth', 1);
title('Permeability (mD)')
axis tight;
view(35, 40);
colorbar('SouthOutside');

%% Initialize schedule and system before solving for all timesteps
%
schedule = deck.SCHEDULE;
system = initADISystem(deck, G, rock, fluid, 'cpr', true);
timer = tic;
[wellSols, states, iter] = runScheduleADI(state, G, rock, system, schedule);
toc(timer)


%% Plot Producer Gas/Oil ratio
%
T = convertTo(cumsum(schedule.step.val), day);
inj  = find([wellSols{1}.sign] == 1);
prod = find([wellSols{1}.sign] == -1);
% Put the well solution data into a format more suitable for plotting
[qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);
figure(1)
clf
gor = qGs(:,prod) ./ qOs(:,prod);
plot(T, gor, '-*b')
xlabel('Days')
title('Gas rate / Oil rate')

%% Plot of Bottom Hole Pressure and gas production/injection rate
% The wells are controlled by gas rate but also constrained by pressure. For
% the producing well, the lower limit for the bottom hole pressure is 1050
% psia (72.395 bar). We observe that the producing well switches control at
% t=5840 days.

figure(2)
clf
subplot(2,2,1)
bhp_p = convertTo(bhp(:,prod), barsa);
plot(T, bhp_p, '-*b')
xlabel('Days')
ylabel('bar')
title('Bottom hole pressure (Producer)')

subplot(2,2,2)
bhp_i = convertTo(bhp(:,inj), barsa);
plot(T, bhp_i, '-*b')
xlabel('Days')
ylabel('bar')
title('Bottom hole pressure (Injector)')

subplot(2,2,3)
qGs_p = convertTo(qGs(:,prod), meter^3/day);
plot(T, qGs_p, '-*b')
xlabel('Days')
ylabel('m^3/day')
title('Gas production rate (m^3/day)')

subplot(2,2,4)
qGs_i = convertTo(qGs(:,inj), meter^3/day);
plot(T, qGs_i, '-*b')
xlabel('Days')
ylabel('m^3/day')
title('Gas injection rate (m^3/day)')

% Compute first time producing well switches to bhp control
ct = cellfun(@(w) w(2).type, wellSols, 'UniformOutput', false);
ind = find(strcmp('bhp', ct), 1);
subplot(2,2,1)
text(T(ind), bhp_p(ind), ...
     'Well control switches from gas rate to bhp \rightarrow  ', ...
     'horizontalalignment', 'right');
subplot(2,2,3)

text(T(ind), qGs_p(ind), ...
     'Well control switches from gas rate to bhp \rightarrow  ', ...
     'horizontalalignment', 'right');

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
