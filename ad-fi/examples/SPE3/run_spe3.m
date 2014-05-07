%% SPE3 case for fully implicit black oil solver
% This example solves the SPE1 problem which consists of gas injection in a
% small ($9\times9\times4$) reservoir with a single producer and
% injector. The problem is parsed and solved from the problem file
% "SPE3.DATA" 

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

%%
% The initial state is a pressure field and  a
%  mixture of water  and free gas (So=0.88) which are constant in each layer. The gas cannot
%  dissolve in oil but ith no initial free
% gas (Sg=0.0) and a constant dissolved gas/oil ratio ("Rs") throughout the
% model.

p0  = deck.SOLUTION.PRESSURE;
sw0 = deck.SOLUTION.SWAT;
sg0 = deck.SOLUTION.SGAS;
s0  = [sw0, 1-sw0-sg0, sg0];
rv0 = deck.SOLUTION.RV;
rs0 = 0;
state = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);   
clear k p0 s0 rv0 rs0

%fluid.rsSat = @(p)rs0;
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

