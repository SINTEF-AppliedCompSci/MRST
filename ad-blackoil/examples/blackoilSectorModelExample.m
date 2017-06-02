%% Example demonstrating use of boundary conditions for pressure support
% We take the SPE1 fluid model to get a simple blackoil-model. We make the
% aqueous phase mobile by manually setting the relative permeability.
[~, ~, fluid, deck, state] = setupSPE1();
fluid.krW = @(s) s.^2;

%% Set up a simple grid with an initial state
% We create a small grid with initial water and oil saturation, as well as
% a Rs (dissolved gas in oil-ratio) of 200. This undersaturated reservoir
% does not have any free gas.
G = cartGrid([50, 50, 1], [1000, 1000, 100]*meter);
G = computeGeometry(G);

rock = makeRock(G, 0.3*darcy, 0.3);

% Define initial composition and pressure
s0 = [0.2, 0.8, 0];
p0 = 300*barsa;
state0 = initResSol(G, p0, s0);
state0.rs = repmat(200, G.cells.num, 1);
% Black oil with disgas
model = ThreePhaseBlackOilModel(G, rock, fluid, 'disgas', true);
%% Set up driving forces
% We will operate a single producer well at a fixed rate. In addition, we
% define a set of boundary conditions at the vertical boundary of the
% domain with a fixed pressure and composition equal to the initial values
% of the reservoir itself.

% Time horizon
T = 10*year;
% Produce 1/4 PV at surface conditions
prate = -0.25*sum(model.operators.pv)/T;
% Oil-rate controlled producer
I = ceil(G.cartDims(1)/2);
J = ceil(G.cartDims(1)/2);
W = [];
W = verticalWell(W, G, rock, I, J, [], ...
    'val',      prate, ...
    'type',     'orat',...
    'comp_i',   [1, 1, 1]/3, ...
    'sign',     -1, ...
    'name',     'Producer');
% Define a lower bhp limit so that we have a lower bound on the reservoir
% pressure during the simulation.
W.lims.bhp = 100*barsa;

% Define boundary conditions
bc = [];
sides = {'xmin', 'xmax', 'ymin', 'ymax'};
for side = 1:numel(sides);
    bc = pside(bc, G, sides{side}, p0, 'sat', s0);
end
% We define initial Rs value for the BC. We can either supply a single
% value per interface or one for all interfaces. In this case, we supply a
% single value for all interfaces since there is no variation in the
% initial conditions.
%
% To define dissolution per face, a matrix of dimensions numel(bc.face)x3x3
% would have to be set.
bc.dissolution = [1, 0,   0;... % Water fractions in phases 
                  0, 1,   0; ...% Oil fractions in phases
                  0, 200, 1];   % Gas fractions in phases
% Simple uniform schedule with initial rampup
dt = rampupTimesteps(T, 30*day);
% Define schedule with well and bc
schedule = simpleSchedule(dt, 'W', W, 'bc', bc);
% Define schedle with only well and no bc
schedule_closed = simpleSchedule(dt, 'W', W, 'bc', []);
%% Simulate case with open boundaries
[ws, states] = simulateScheduleAD(state0, model, schedule);

%% Simulate case with closed boundaries
[ws_c, states_c] = simulateScheduleAD(state0, model, schedule_closed);

%% Compare the two results
% We observe that the simulation with open boundary has a very different
% pressure build up than the same problem with closed boundaries. Note also
% that once the bhp limit is reached, the well switches controls in the
% closed model, and the oil production rate changes.
plotWellSols({ws, ws_c}, cumsum(schedule.step.val), ...
    'datasetnames', {'Open boundaries', 'Closed boundaries'})
%% Plot the gas in the well cells
% If we do not have open boundaries, the pressure will eventually drop as a
% result of the large removed volumes of oil and dissolved gas. For the
% model with closed boundaries, we get drop-out of free gas.
gas_open = cellfun(@(x) x.s(W.cells(1), 3), states);
gas_closed = cellfun(@(x) x.s(W.cells(1), 3), states_c);

clf;
plot([gas_open, gas_closed], 'linewidth', 2)
ylabel('S_g')
xlabel('Step #');
legend('Open boundaries', 'Closed boundaries');