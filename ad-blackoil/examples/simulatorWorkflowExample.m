%% Workflow example for MRST-AD
% This example aims to show the complete workflow for creating, running and
% analyzing a simulation model. Unlike the other examples, we will create
% all features of the model manually to get a self-contained script without
% any input files required.
%
% The model we will setup is a two phase oil/water model with some
% compressibility and multiple wells. We begin by setting up the grid and
% rock structures. The grid is created by "makeModel3", which creates a
% structured model with intersecting faults.
%
% We also define a layered permeability structure, with 300, 100 and 500
% darcy in the lower, middle and top layers respectively.
%
% Note that this example shows a simple conceptual model designed to show
% the workflow rather than a problem representing a realistic scenario in
% terms of well locations and fluid physics.

% Define grid
grdecl = makeModel3([50, 50, 5], [1000, 1000, 5]*meter);
G = processGRDECL(grdecl);
G = computeGeometry(G);

% Set up permeability based on K-indices
[I, J, K] = gridLogicalIndices(G);

top = K < G.cartDims(3)/3;
lower = K > 2*G.cartDims(3)/3;
middle = ~(lower | top);

px = ones(G.cells.num, 1);
px(lower) = 300*milli*darcy;
px(middle) = 100*milli*darcy;
px(top) = 500*milli*darcy;

% Introduce anisotropy by setting K_x = 10*K_z.
perm = [px, px, 0.1*px];
rock = makeRock(G, perm, 0.3);
%% Define wells and simulation schedule
% We will insert three producers, operating at fixed bottom hole pressure
% and perforated throughout all layers of the model along with a single
% injector injecting one pore volume over 10 years (the total simulation
% length). We also set up a schedule consisting of 5 small control steps
% initially, followed by 25 larger steps. We keep the well controls fixed
% throughout the simulation.

simTime = 10*year;
nstep = 25;
refine = 5;

pv = poreVolume(G, rock);
injRate = 1*sum(pv)/simTime;

% Place wells
offset = 10;
W = [];
% Producers
W = verticalWell(W, G, rock, offset, offset, [],...
                'Name', 'P1', 'comp_i', [1 0], 'Val', 250*barsa, 'Type', 'bhp');
W = verticalWell(W, G, rock, G.cartDims(1) - offset, G.cartDims(2) - offset, [],...
                'Name', 'P2', 'comp_i', [1 0], 'Val', 250*barsa, 'Type', 'bhp');
W = verticalWell(W, G, rock, offset, G.cartDims(2) - offset, [], ...
                'Name', 'P3', 'comp_i', [1 0], 'Val', 250*barsa, 'Type', 'bhp');
% Injectors
W = verticalWell(W, G, rock, floor(G.cartDims(1)/2) + 4, floor(G.cartDims(1)/2) - 3, 1,...
                'Name', 'Injector', 'comp_i', [1 0], 'Val', injRate, 'Type', 'rate');

% Compute the timesteps
startSteps = repmat((simTime/(nstep + 1))/refine, refine, 1);
restSteps =  repmat(simTime/(nstep + 1), nstep, 1);
timesteps = [startSteps; restSteps];
% Set up the schedule containing both the wells and the timesteps
schedule = simpleSchedule(timesteps, 'Wells', W);

% Plot horizontal permeability and wells
clf;
plotCellData(G, rock.perm(:, 1)/(milli*darcy))
plotWell(G, W)
axis tight
view(50, 50)
colorbar
title('K_x (mD)')
%% Set up simulation model
% We set up a two phase simulation model based on automatic
% differentiation. We define a three phase fluid (with properties given for
% oil, water and gas) but the model will only use the oil and water values.
% Once we have defined the fluid, we further modify the oil compressibility
% by changing the bO-factor to have some constant compressibility using
% an anonymous function of pressure (using the standard Matlab builtin exp)
% to get a function on the form
%
% $ b_o(p) = b_0 e^{(p - p_0) c} $
%
% We then pass the fundamental structures (grid, rock and fluid) onto the
% two-phase oil/water model constructor.
mrstModule add ad-core ad-blackoil ad-props

fluid = initSimpleADIFluid('mu',    [1, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);

c = 0.001/barsa;
p_ref = 300*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);

clf;
p0 = (100:10:500)*barsa;
plot(p0/barsa, fluid.bO(p0))
xlabel('Pressure (bar)')
ylabel('Reference to reservoir density ratio for the oil phase (bO)')
gravity reset on
model = TwoPhaseOilWaterModel(G, rock, fluid);
%% Define initial state
% Once we have a model, we need to set up a initial state. We set up a very
% simple initial state. MRST uses water, oil, gas ordering internally, so
% in this case we have water in the first column and oil in the second for
% the saturations. We let the bottom part of the reservoir be completely
% water filled, and the top completely oil filled.
sW = ones(G.cells.num, 1);
sW(G.cells.centroids(:, 3) < 5) = 0;

sat = [sW, 1 - sW];

g = model.gravity(3);
% Compute initial pressure
p_res = p_ref + g*G.cells.centroids(:, 3).*[sW.*model.fluid.rhoWS + (1 - sW).*model.fluid.rhoOS];
state0 = initResSol(G, p_res, sat);
%% Simulate base case
% Once we have defined the schedule (dynamic controls and time), model
% (mathematical description of how to advance the solution) and the initial
% solution, we can simulate the problem.
mrstVerbose on
[wellSols, states] = simulateScheduleAD(state0, model, schedule);
%% Plot results
% We launch a plotting tool for both the reservoir quantities (pressures
% and saturations, located in states) and the well solution (well rates and
% bottom hole pressures, located in wellSols).
mrstModule add mrst-gui
clf
plotToolbar(G, states)
view(50, 50);
% Wells
simtime = cumsum(schedule.step.val);
plotWellSols(wellSols, simtime, 'field', 'qOs');

%% Create a upscaled, coarser model
% The fine scale model has approximately 10000 cells. If we want a smaller
% model we can easily define an upscaled model. Here we set up a simple
% uniform partition of approximately 50 cells based on the logical
% partition.
mrstModule add coarsegrid
cdims = [5, 5, 2];
p0 = partitionUI(G, cdims);

clf
plotCellData(G, mod(p0, 13), 'EdgeColor', 'none')
axis tight off
view(50, 50)
%% Split blocks over the faultlines
% We see that several coarse blocks cross the fault lines. To get
% hexahedral coarse blocks, we create a grid where the faults act as
% barriers and apply the "processPartition" routine to split any coarse
% blocks intersected by faults.
%
% Afterwards, we show the new partition and highlight blocks created due to
% the modification of the fault.

G_fault = makeInternalBoundary(G, find(G.faces.tag > 0));
p = processPartition(G_fault, p0);

clf
plotCellData(G, mod(p, 13), 'EdgeColor', 'none')
plotGrid(G, p ~= p0, 'EdgeColor', 'r', 'LineWidth', 2, 'FaceColor', 'none')
axis tight off
view(50, 50)
title('Partition after processing by faults');

%% Upscale the model and run the coarser problem
% We can now directly upscale the model, schedule and initial state
% automatically. By default, the upscaling routines uses the simplest
% possible options, i.e. harmonic averaging of permeabilities. It is
% possible to use more advanced options, but for the purpose of this
% example we will use the defaults.
%
% Once we have a upscaled model, we can again simulate the new schedule and
% observe that the time taken is greatly reduced.
model_c = upscaleModelTPFA(model, p);
G_c    = model_c.G;
rock_c = model_c.rock;

schedule_c = upscaleSchedule(model_c, schedule);
state0_c = upscaleState(model_c, model, state0);

[wellSols_c, states_c] = simulateScheduleAD(state0_c, model_c, schedule_c);

%% Plot the coarse results, and compare the well solutions
% We plot the coarse solutions and compare the well solutions. Note that
% the upscaling will result in only 70 cells, which is unlikely to give
% good results with only simple harmonic averaging of permeabilities.
clf
plotToolbar(G_c, states_c)
view(50, 50);

plotWellSols({wellSols, wellSols_c}, simtime, 'DatasetNames', {'Fine scale', 'Upscaled'}, 'Field', 'qOs');
%% Compute flow diagnostics
% As an alternative to looking at well curves, we can also look at the flow
% diagnostics of the models. Flow diagnostics are simple routines based on
% time-of-flight and tracer equations, which aim to give a qualitative
% understanding of reservoir dynamics. Here, we take the end-of-simulation
% states as a snapshot for both the fine and coarse model and compute
% time-of-flight and well tracers.
mrstModule add diagnostics
D   = computeTOFandTracer(states{end},   G,   rock,   'Wells', schedule.control.W);
D_c = computeTOFandTracer(states_c{end}, G_c, rock_c, 'Wells', schedule_c.control.W);
%% Plot total arrival times
% We plot the sum of forward and backwards time of flight, showing which
% regions have large amounts of flow. Large values indicate that a cell
% has little flow going through it.
%
% Since the values vary by several orders of magnitude, we take the
% $log_{10}$ transform before plotting. We also use the same color axis to
% ensure that the plots can be compared.
figure(1); clf
plotCellData(G, log10(sum(D.tof, 2)));
view(50, 50);
title('Log of total travel time, fine model');
c = caxis();

figure(2); clf
plotCellData(G_c, log10(sum(D_c.tof, 2)));
view(50, 50);
title('Log of total travel time, coarse model');
caxis(c)
%% Plot tracer partitioning
% We can also look at the tracer partitioning for the producers, showing
% the drainage regions for the different wells.
%
% See the diagnostics module for more examples and more in-depth
% discussions of how flow diagnostics can be used.
figure(1); clf
plotCellData(G, D.ppart);
view(50, 50);
title('Drainage regions, fine model');

figure(2); clf
plotCellData(G_c, D_c.ppart);
view(50, 50);
title('Drainage regions, coarse model');
%% Launch interactive diagnostics tools
% We can also examine the diagnostics interactively using the diagnostics
% viewer.
close all;
interactiveDiagnostics(G, rock, schedule.control.W, 'state', states{end}, 'computeFlux', false);
interactiveDiagnostics(G_c, rock_c, schedule_c.control.W, 'state', states_c{end}, 'computeFlux', false);