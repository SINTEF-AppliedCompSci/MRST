%% Classic Buckley-Leverett Problem: 1D Horizontal Displacement
% This example uses the classical Buckley-Leverett problem to introduce you
% to basic functionality in the object-oriented automatic differentiation
% framework.
%
% The Buckley-Leverett problem describes an incompressible displacement in
% a 1D homogeneous and horizontal medium with water injected at the left
% boundary x=0 and fluids produced at the right boundary x=L. We assume
% that the fluids have unit viscosity and density and have relative
% permeabilities that follow a standard quadratic Corey model

mrstModule add ad-blackoil ad-core ad-props mrst-gui

%% Set up model
% We start by generating a model object that describes the reservoir. To
% construct this model, we need three more fundamental structures: 'G'
% represents the grid with reservoir geometry, 'rock' holds the
% petrophysical properties, and 'fluid' gives the fluid properties. In
% addition to the model object, we need a 'state' struct that defines the
% initial state in the reservoir (pressure and fluid saturations and
% compositions).

% Construct 3D grid with 50 cells in the x-direction
G = cartGrid([50, 1, 1], [1000, 10, 10]*meter);
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', darcy*ones(G.cells.num, 1), ...
              'poro', .3*ones(G.cells.num, 1));

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'WO', 'n', [2 2]);

% Set up model and initial state.
model = TwoPhaseOilWaterModel(G, rock, fluid);
state0 = initResSol(G, 50*barsa, [0, 1]);
state0.wellSol = initWellSolAD([], model, state0);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
pv = poreVolume(G, rock);
injRate = -sum(pv)/(500*day);
bc = fluxside([], G, 'xmin', -injRate, 'sat', [1, 0]);
bc = pside(bc, G, 'xmax', 0*barsa, 'sat', [0, 1]);

%% Simulate 1 PVI using a manual loop
% There are several ways to run a simulation. A simple approach is to use a
% manual loop, where you explicitly call a nonlinear solver to solve each
% implicit time step
solver = NonLinearSolver();

dT = 20*day;
n = 25;
states = cell(n+1, 1);
states{1} = state0;
for i = 1:n
    state = solver.solveTimestep(states{i}, dT, model, 'bc', bc);
    states{i+1} = state;
end

%% Plot the result
% We set up a plot using plotToolbar from the mrst-gui module. Since the
% problem is essentially one dimensional, we plot the water saturation
% (first column of the "s" field in state) as a 1D plot.
close all
plotToolbar(G, states, 'field', 's:1', 'plot1d', true, ...
                       'lockCaxis', true, 'startplayback', true);

%% Repeat simulation with general solver
% To use the general solver, we first need to set up a schedule that
% describes the time steps and the drive mechanisms (wells, boundary
% conditions, and source terms) that are active in each time step. In
% addition, one can specify various forms of time-step control. Here,
% however, we simply rely on the default setup
schedule = simpleSchedule(repmat(dT,1,25), 'bc', bc);
[~,sstates] = simulateScheduleAD(state0, model, schedule);

close all
plotToolbar(G, sstates, 'field', 's:1','lockCaxis',true), 
caxis([0 1]), view(10,10)
colorbar

%% Repeat simulation with visualization
% The general solver has a hook, that enables you to visualize the progress
% of the simulation (and stop it and continue running it in 'debug' mode).
close all    
fn = getPlotAfterStep(state0, model, schedule, ...
    'plotWell', false, 'plotReservoir', true, 'field', 's:1', ...
    'lockCaxis',true, 'plot1d', true);

[~,sstates,report] = ...
   simulateScheduleAD(state0, model, schedule,'afterStepFn', fn);
