%% Set up case
%
% We set up a very simple incompressible model with boundary conditions
mrstModule add ad-blackoil ad-core ad-props mrst-gui

G = cartGrid([10, 1, 1], [100, 10, 10]*meter);
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', 1000*ones(G.cells.num, 1)*milli*darcy, ...
              'poro', ones(G.cells.num, 1));

% Default fluid with unit values
fluid = initSimpleADIFluid();

% Set up model and initial state.
model = TwoPhaseOilWaterModel(G, rock, fluid);
state0 = initResSol(G, 50*barsa, [1, 0]);
state0.wellSol = initWellSolAD([], model, state0);

bc = [];
src = [];

pv = poreVolume(G, rock);
injRate = -sum(pv)/(500*day);

bc = fluxside(bc, G, 'xmin', -injRate, 'sat', [.25, .75]);
bc = pside(bc, G, 'xmax', 0*barsa, 'sat', [.25, .75]);

%% Simulate some timesteps
solver = NonLinearSolver();

dT = 50*day;
n = 20;
states = cell(n+1, 1);
states{1} = state0;
for i = 1:n
    state = solver.solveTimestep(states{i}, dT, model, 'bc', bc, 'src', src);
    states{i+1} = state;
end

%% Plot the result
close all
plotToolbar(G, states)
colorbar