%% SPE3 case using fully implicit black oil solver
% This example solves the SPE3 problem which consists of gas injection in a
% small ($9\times9\times4$) reservoir. The problem is originally a compositional
% problem. Using PVTi, we convert it to a blackoil problem and the resulting
% datas are provided in the file "SPE3.DATA". The oil can vaporize but the gas
% cannot dissolve in oil so that the gas/oil ratio remains equal to zero during
% the whole simulation. 

require ad-fi deckformat

% Read and process file.
current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'SPE3.DATA');

deck = readEclipseDeck(fn);

% The deck is given in field units, MRST uses metric.
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

% Create a special ADI fluid which can produce differentiated fluid
% properties.
fluid = initDeckADIFluid(deck);

% The case includes gravity
gravity on

%% Set up initial state
%  The initial state is a pressure field and a mixture of water and free gas
%  corresponding to the equilibrium between gravitational and capillary forces.

p0  = deck.SOLUTION.PRESSURE;
sw0 = deck.SOLUTION.SWAT;
sg0 = deck.SOLUTION.SGAS;
s0  = [sw0, 1-sw0-sg0, sg0];
rv0 = deck.SOLUTION.RV;
rs0 = 0;
state = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);   
clear k p0 s0 rv0 rs0

%% Plot well and permeability
% The permeability is constant in each layer. There is one injecting and one
% producing well.

clf;
W = processWells(G, rock, deck.SCHEDULE.control(1));
plotCellData(G, convertTo(rock.perm(:,1), milli*darcy), 'FaceAlpha', .5, ...
            'EdgeAlpha', .3, 'EdgeColor', 'k');
plotWell(G, W, 'fontsize', 10, 'linewidth', 1);
title('Permeability (mD)')
axis tight;
view(35, 40);
colorbar('SouthOutside');

%% Initialize schedule and system before solving for all timesteps
schedule = deck.SCHEDULE;
system = initADISystem(deck, G, rock, fluid, 'cpr', true);
timer = tic;
[wellSols, states, iter] = runScheduleADI(state, G, rock, system, schedule);
toc(timer)  

%% Plot Producer Gas/Oil ratio
T = convertTo(cumsum(schedule.step.val), day);
inj  = find([wellSols{1}.sign] == 1);
prod = find([wellSols{1}.sign] == -1);
% Put the well solution data into a format more suitable for plotting
[qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);
figure(1)
clf
gor = qGs(:,prod)./qOs(:,prod);
plot(T, gor, '-*b')
xlabel('Days')
title('Gas rate / Oil rate')


%% Plot Producer Bottom Hole Pressure
figure(2)
clf
bhp_p = bhp(:,prod);
plot(T,     convertTo(bhp_p, barsa), '-*b')
xlabel('Days')
ylabel('bar')
title('Bottom hole pressure (Producer)')


%% Plot Injector Bottom Hole Pressure
figure(3)
clf
bhp_i = bhp(:,inj);
plot(T,     convertTo(bhp_i, barsa), '-*b')
xlabel('Days')
ylabel('bar')
title('Bottom hole pressure (Injector)')

