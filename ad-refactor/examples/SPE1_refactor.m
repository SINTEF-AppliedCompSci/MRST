%% SPE1 case for fully implicit black oil solver
% This example solves the SPE1 problem which consists of gas injection in a
% small ($10\times10\times3$) reservoir with a single producer and
% injector. The problem is parsed and solved from the problem file
% "odeh_adi" and the result is then compared to output from a major
% commercial reservoir simulator (Eclipse 100).

try
   require ad-fi deckformat mrst-gui ad-blackoil
catch
   mrstModule add ad-fi deckformat mrst-gui ad-blackoil
end

% Read and process file.
current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'odeh_adi.data');

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
gravity reset on


% The initial state is a pressure field that is constant in each layer, a
% uniform mixture of water (Sw=0.12) and oil (So=0.88) with no initial free
% gas (Sg=0.0) and a constant dissolved gas/oil ratio ("Rs") throughout the
% model.
%
% The pressure and Rs values are derived through external means.
clear prod
[k, k] = ind2sub([prod(G.cartDims(1:2)), G.cartDims(3)], ...
                  G.cells.indexMap);  %#ok

p0    = [ 329.7832774859256 ; ...  % Top layer
          330.2313357125603 ; ...  % Middle layer
          330.9483500720813 ];     % Bottom layer

p0    = convertFrom(p0(k), barsa);
s0    = repmat([ 0.12, 0.88, 0.0 ], [G.cells.num, 1]);
rs0   = repmat( 226.1966570852417 , [G.cells.num, 1]);
rv0   = 0; % dry gas

state = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);   clear k p0 s0 rs0

% Determine the model automatically from the deck. It should be a
% three-phase black oil model with gas dissoluton.
model = selectModelFromDeck(G, rock, fluid, deck);

model %#ok, intentional display

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(G, model, rock, deck);

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

plotCellData(G, convertTo(rock.perm(:,1), milli*darcy), 'FaceAlpha', .5, ...
            'EdgeAlpha', .3, 'EdgeColor', 'k');
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
xlabel('Years')
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
xlabel('Years')
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
xlabel('Years')
ylabel('bar')
title('Bottom hole pressure (Injector)')
