%% Example demonstrating a three dimensional, six component problem
% We set up a simple grid and rock structure. The numbers can be adjusted
% to get a bigger/smaller problem. The default is a small problem with
% 20x20x2 grid blocks.
mrstModule add ad-core compositional ad-props mrst-gui
% Dimensions
nx = 20;
ny = nx;
nz = 2;
% Name of problem and pressure range
casename = 'lumped_1';
minP = 100*barsa;
maxP = 200*barsa;
% Set up grid and rock
dims = [nx, ny, nz];
pdims = [1000, 1000, 1];
G = cartGrid(dims, pdims);
G = computeGeometry(G);
rock = makeRock(G, 50*milli*darcy, 0.25);
%% Set up quarter five spot well pattern
% We place vertical wells in opposing corners of the reservoir. The
% injector is rate controlled and the producer is bottom hole pressure
% controlled.
W = [];
% Injector
W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', [1, 0], 'name', 'Inj',...
    'Val', 0.0015, 'sign', 1, 'type', 'rate');
% Producer
W = verticalWell(W, G, rock, nx, ny, [], ...
    'comp_i', [0.5, 0.5], 'Name', 'Prod', 'Val', minP, 'sign', -1, 'Type', 'bhp');
%% Set up model and initial state
% We set up a problem with quadratic relative permeabilities. The fluid
% model is retrieved from "High Order Upwind Schemes for Two-Phase,
% Multicomponent Flow" (SPE 79691) by B. T. Mallison et al.
%
% The model consists of six components. Several of the components are not
% distinct molecules, but rather lumped groups of hydrocarbons with similar
% molecular weight. The reservoir contains all these components initially,
% which is then displaced by the injection of pure CO2.

nkr = 2;
[fluid, info] = getCompositionalFluidCase(casename);
flowfluid = initSimpleADIFluid('n', [nkr, nkr, nkr], 'rho', [1000, 800, 10]);

gravity reset on
model = ThreePhaseCompositionalModel(G, rock, flowfluid, fluid, 'water', false);
% Reduce tolerances somewhat, as the defaults are very strict.
model.incTolPressure = 1e-2;
model.toleranceWellRate = 1e-2;

ncomp = fluid.getNumberOfComponents();
state0 = initResSol(G, minP + (maxP - minP)/2, 1);

state0.T = repmat(info.temp, G.cells.num, 1);
state0.components = cell(1, ncomp);
for i = 1:ncomp
   state0.components{i} = repmat(info.initial(i), G.cells.num, 1);
end
state0.wellSol = initWellSolAD(W, model, state0);
state0 = model.computeFlash(state0, inf);

for i = 1:numel(W)
    W(i).components = info.injection;
end
%% Set up schedule and simulate the problem
% We simulate two years of production with a geometric rampup in the
% timesteps.
time = 2*year;
n = 45;
dt = rampupTimesteps(time, 7*day, 5);
schedule = simpleSchedule(dt, 'W', W);

[ws, states, rep] = simulateScheduleAD(state0, model, schedule);
%% Plot all the results
lf = get(0, 'DefaultFigurePosition');
h = figure('Position', lf + [0, -200, 350, 200]);
nm = ceil(ncomp/2);
v = [-30, 60];
for step = 1:numel(states)
    figure(h); clf
    state = states{step};
    for i = 1:ncomp
        subplot(nm, 3, i);
        plotCellData(G, state.components{i}, 'EdgeColor', 'none');
        view(v);
        title(fluid.names{i})
        caxis([0, 1])
    end
    subplot(nm, 3, ncomp + 1);
    plotCellData(G, state.pressure, 'EdgeColor', 'none');
    view(v);
    title('Pressure')
    
    subplot(nm, 3, ncomp + 2);
    plotCellData(G, state.s(:, 1), 'EdgeColor', 'none');
    view(v);
    title('sO')
    
    subplot(nm, 3, ncomp + 3);
    plotCellData(G, state.s(:, 2), 'EdgeColor', 'none');
    view(v);
    title('sG')
    drawnow
end
%% Plot the results in the interactive viewer
figure(1); clf;
plotToolbar(G, states)
view(v);
axis tight