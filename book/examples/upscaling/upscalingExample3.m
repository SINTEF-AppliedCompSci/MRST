%% Upscaling workflow
% This example aims to show the complete workflow for creating, running and
% analyzing a simulation model. Unlike the other examples, we will create
% all features of the model manually to get a self-contained script without
% any input files required.
%
% The model we setup is a slightly compressible two-phase oil/water model
% with multiple wells. The reservoir has a layered strategraphy consisting
% of three layered sections, four layers with different permeability, and
% three major faults.
%
% Note that this example features a simple conceptual model designed to
% show the workflow rather than a problem representing a realistic scenario
% in terms of well locations and fluid physics.
mrstModule add ad-core ad-blackoil ad-props diagnostics mrst-gui
close all;

%% Horizons for internal and external geology
% We begin by building the horizons that define the top and bottom
% structure of the sector, as well as two internal erosions.

% Define areal mesh
[xmax,ymax, n]  = deal(1000*meter, 1000*meter, 30);
[x, y] = meshgrid(linspace(0,xmax,n+1), linspace(0,ymax,n+1));
[x, y] = deal(x',y');

% Basic dome structure
dome = 1-exp(sqrt((x - xmax/2).^2 + (y - ymax/2).^2)*1e-3);

% Periodic perturbation
[xn,yn] = deal(pi*x/xmax,pi*y/ymax);
perturb = sin(5*xn) + .5*sin(4*xn+6*yn) + cos(.25*xn*yn./pi^2) + cos(3*yn);
perturb = perturb/3.5;

% Random small-scale perturbation
rng(0);
[h, hr] = deal(8,1);
zt = 50 + h*perturb + rand(size(x))*hr - 20*dome;
zb = zt + 30;
zmb = min(zb + 4 + 0.01*x - 0.020*y + hr*rand(size(x)), zb);
zmt = max(zb -15 + 0.01*x - 0.025*y + hr*rand(size(x)), zt);
horizons = {struct('x', x, 'y', y, 'z', zt), ...
    struct('x', x, 'y', y, 'z', zmt), ...
    struct('x', x, 'y', y, 'z', zmb), ...
    struct('x', x, 'y', y, 'z', zb)};

surf(x,y,zt-.2, 'EdgeC','r','FaceC',[.8 .8 .8]),  hold on
mesh(x,y,zmt-.2,'EdgeC','g','FaceC',[.7 .7 .7]), 
mesh(x,y,zmb-.2,'EdgeC','g','FaceC',[.6 .6 .6]), 
mesh(x,y,zb-.2, 'EdgeC','b','FaceC',[.5 .5 .5]); hold off
set(gca,'ZDir','reverse')
view(-50,10); axis off

%% Interpolate to build unfaulted corner-point grid
dims = [40, 40]; layers = [3 6 3];
grdecl = convertHorizonsToGrid(horizons, 'dims', dims, 'layers', layers);
G = processGRDECL(grdecl);
[~,~,k] = gridLogicalIndices(G);

figure, plotCellData(G,k,'EdgeAlpha',.2); view(3);
colormap(.5*jet + .5*ones(size(jet)));
view(-50,30); axis off

%% Insert faults
[X,Y,Z]  = buildCornerPtNodes(grdecl);

i=47:80; Z(i,:,:) = Z(i,:,:) + .022*min(0,Y(i,:,:)-550);
j= 1:30; Z(:,j,:) = Z(:,j,:) + .021*min(0,X(:,j,:)-400);
j=57:80; Z(:,j,:) = Z(:,j,:) + .023*min(0,X(:,j,:)-750);
grdecl.ZCORN = Z(:);

G = processGRDECL(grdecl);
G = computeGeometry(G);
[~,~,k] = gridLogicalIndices(G);

figure, plotCellData(G,k,'EdgeAlpha',.2); view(3);
plotFaces(G,find(G.faces.tag>0),'EdgeColor','r','FaceColor',[.8 .8 .8]);
colormap(.5*jet + .5*ones(size(jet)));
view(-50,30); axis tight off

%% Petrophysics
% Set up permeability based on K-indices and introduce anisotropy by
% setting K_z = .1*K_x
rng(357371);
[K,L] = logNormLayers(G.cartDims, [100 400 10 50]*milli*darcy);
K = K(G.cells.indexMap);
perm = [K, K, 0.1*K];
rock = makeRock(G, perm, 0.3);

% Plot horizontal permeability
figure
K = convertTo(K,milli*darcy);
plotCellData(G, log10(K),'EdgeAlpha',.1)
colorbarHist(K,[1 1500],'South',100,true);
view(-50, 50), axis tight off

%% Define wells
% Hydrocarbon is recovered from producers, operating at fixed bottom-hole
% pressure and perforated throughout all layers of the model. The producers
% are supported by a single water injector set to inject one pore volume
% over 10 years (the total simulation length).

% Producers
simTime = 10*year;
pv      = poreVolume(G, rock);
injRate = 1*sum(pv)/simTime;
offset  = 10;
W = verticalWell([], G, rock, offset, offset, [],...
                'Name', 'P1', 'comp_i', [1 0], ...
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50);
W = verticalWell(W, G, rock,  offset, floor(G.cartDims(1)/2)+3, [],...
                'Name', 'P2', 'comp_i', [1 0], ...
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50);
W = verticalWell(W, G, rock, offset, G.cartDims(2) - offset/2, [], ...
                'Name', 'P3', 'comp_i', [1 0], ...
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50);

% Injectors
W = verticalWell(W, G, rock, G.cartDims(1)-5, offset, [],...
                'Name', 'I1', 'comp_i', [1 0], ...
                'Val', injRate, 'Type', 'rate', 'refDepth', 50);

% Plot the wells
plotWell(G, W,'color','k')
axis tight

%% Define fluid behavior and instantiate model
% We set up a two-phase oil-water simulation model based on automatic
% differentiation. The resulting object is a special case of a general
% three-phase model and to instantiate it, we start by constructing a
% three-phase fluid structure with properties given for oil, water, and
% gas. (The gas properties will be neglected once we construct the
% two-phase object). Water is assumed to be incompressible, whereas oil has
% constant compressibility, giving an expansion factor of the form, $b_o(p)
% = b_0 exp[c (p - p_0)]$. To define this relation, we set the 'bo' field
% of the fluid structure to be an anonymous function that calls the builtin
% 'exp' function with appropriate arguments. Since the fluid model is a
% struct containing function handles, it is simple to modify the fluid to
% use alternate functions. We then pass the fundamental structures (grid,
% rock and fluid) onto the two-phase oil/water model constructor.

% Three-phase template model
fluid = initSimpleADIFluid('mu',    [1, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);

% Constant oil compressibility
c        = 0.001/barsa;
p_ref    = 300*barsa;
fluid.bO = @(p, varargin) exp((p - p_ref)*c);

% Construct reservoir model
gravity reset on
model = TwoPhaseOilWaterModel(G, rock, fluid);


%% Define initial state
% Once we have a model, we need to set up an initial state. We set up a
% very simple initial state; we let the bottom part of the reservoir be
% completely water filled, and the top completely oil filled. MRST uses
% water, oil, gas ordering internally, so in this case we have water in the
% first column and oil in the second for the saturations.

depthOW = 85*meter;
depthD  = 10*meter;
region = getInitializationRegionsBlackOil(model, depthOW, ...
            'datum_depth', depthD, 'datum_pressure', p_ref);
state0 = initStateBlackOilAD(model, region);

figure
plotCellData(G, state0.s(:,1),'EdgeAlpha',.1), colormap(flipud(winter));
plotWell(G,W,'color','k')
patch([-50 1050 1050 -50],[-50 -50 1050 1050],depthOW*ones(1,4), ...
    ones(1,4), 'FaceColor',[.6 .6 1],'EdgeColor','r','LineWidth',1);
view(-50, 50), axis tight off

%% Define simulation schedule and set solver parameters
% We define a relatively simple schedule consisting of five small control
% steps initially, followed by 25 larger steps. We keep the well controls
% fixed throughout the simulation. To accelerate the simulation, we set
% somewhat stricter tolerances and use a CPR preconditioner with an
% algebraic multigrid solver (AMGCL) for the elliptic pressure system

% Compute the timestep
nstep   = 25;
refine  = 5;
startSteps = repmat((simTime/(nstep + 1))/refine, refine, 1);
restSteps=  repmat(simTime/(nstep + 1), nstep, 1);
timesteps = [startSteps; restSteps];

% Set up the schedule containing both the wells and the timestep
schedule = simpleSchedule(timesteps, 'W', W);

% Tighten tolerences
model.drsMaxRel = inf;
model.dpMaxRel  = .1;
model.dsMaxAbs  = .1;

% Set up CPR preconditioner
try
    mrstModule add linearsolvers
    callAMGCL(speye(5), ones(5,1));
    pressureSolver = AMGCLSolverAD('tolerance', 1e-4);
catch
    warning(['Failed to use algebraic multigrid solver from AMGCL. ' ...
        'Trying to proceed with the default direct solver in MATLAB.']);
    pressureSolver = BackslashSolverAD();
end
linsolve = CPRSolverAD('ellipticSolver', pressureSolver, 'relativeTolerance', 1e-3);


%% Simulate base case
% We now have an intial state, the schedule defines dynamic controls and
% time steps, and the model gives the mathematical description of how to
% advance the solution one time step. We then have all we need to simulate
% the problem. Because the simulation will consume some time, we launch a
% progress report and a plotting tool for the well solutions (well rates
% and bottom-hole pressures)
fn = getPlotAfterStep(state0, model, schedule,'view',[50 50], ...
                     'field','s:1','wells',W);
[wellSols, states, report] = ...
   simulateScheduleAD(state0, model, schedule, ...
                       'LinearSolver', linsolve, 'afterStepFn',fn);


%% Plot reservoir states
% We launch a plotting tool for the reservoir quantities (pressures
% and saturations, located in states).
figure, plotToolbar(G, states)
view(50, 50);
plotWell(G,W);

figure
mrstModule add coarsegrid
CG = generateCoarseGrid(G,ones(G.cells.num,1));
plotWell(G, W,'color','k')
plotFaces(CG,1:CG.faces.num,'EdgeColor','k','FaceColor','none','LineWidth',2);
s=states{6}.s(:,1); plotCellData(G,s,s>1e-2); caxis([0 1]);
colormap(flipud(winter)), view(-53,46); axis tight off

%% Launch flow diagnostics
interactiveDiagnostics(G, rock, W, 'state', states{end}, 'computeFlux', false);
set(gcf,'position',[440 317 866 480]); axis normal
view(-85,73)

%% Create an upscaled, coarser model
% The fine scale model has 16 350 cells. If we want a smaller model we can
% easily define an upscaled model. We start by setting up a straightforward
% uniform partition of the IJK-indices. This gives 1388 coarse blocks
% (after a postprocessing to split disconnected blocks), which corresponds
% to more than 20 times reduction in the number of unknowns.
mrstModule add coarsegrid
cdims = [10, 10, 12];
p0 = partitionUI(G, cdims);
% p0 = partitionLayers(G, cdims(1:2), L);

figure;
CG = generateCoarseGrid(G,p0);
plotGrid(G,'FaceColor',[1 .8 .9],'EdgeAlpha',.1);
plotFaces(CG,1:CG.faces.num,'EdgeColor','k','FaceColor','none','LineWidth',1.5);
axis tight off, view(125, 55)
title('Straightforward index partition');

%% Split blocks over the fault lines
% The partition above may give blocks that have cells on opposite sides of
% a fault. To ensure that as many as possible of the coarse blocks are
% hexahedral (except for those that are partially eroded), we split blocks
% that have cells on both sides of a fault. To do this, we introduce a
% temporary grid in which the faults act as barriers, and then perform a
% simple postprocessing to split any coarse blocks intersected by one of
% the faults.
%
% Afterwards, we show the new partition and highlight blocks created due to
% the modification of the fault.
Gf = makeInternalBoundary(G, find(G.faces.tag > 0));
p = processPartition(Gf, p0);
plotGrid(G,p>prod(max(p0)),'FaceColor','g','EdgeColor','none');
title('Splitting over fault lines');

%% Upscale the model and simulate the coarser problem
% We can now directly upscale the model, schedule, and initial state. By
% default, the upscaling routine uses the simplest possible options, i.e.,
% harmonic averaging of permeabilities. It is possible to use more advanced
% options, but for the purpose of this example we will use the defaults.
%
modelC1    = upscaleModelTPFA(model, p);
scheduleC1 = upscaleSchedule(modelC1, schedule);
state0C1   = upscaleState(modelC1, model, state0);

% Plot the initial state
figure,
GC    = modelC1.G;
rockC = modelC1.rock;
plotCellData(GC,state0C1.s(:,1),'EdgeColor','none'), colormap(flipud(winter));
plotWell(G,W,'color','k')
plotFaces(GC,1:GC.faces.num,'FaceColor','none');
patch([-50 1050 1050 -50],[-50 -50 1050 1050],85*ones(1,4), ...
    ones(1,4), 'FaceColor',[.6 .6 1],'EdgeColor','r','LineWidth',1);
view(-50, 50), axis tight off

% Once we have an upscaled model, we can again simulate the new schedule
% and observe that the time taken is greatly reduced, even without the use
% of a CPR preconditioner.
[wellSolsC1, statesC1] = simulateScheduleAD(state0C1, modelC1, scheduleC1);

%% Create a second upscaled model with flow-based upscaling
% We use the upscaling module to create a tailored upscaled model. This
% upscaling routine uses an incompressible flow field with the wells of the
% problem to perform a global upscaling.
mrstModule add incomp agglom upscaling

[~, TC, WC] = upscaleTrans(GC, model.operators.T_all, ...
    'Wells', W, 'bc_method', 'wells', 'fix_trans', true);

modelC2 = upscaleModelTPFA(model, p, 'transCoarse', TC);
scheduleC2 = schedule;
for i = 1:numel(scheduleC2.control)
    scheduleC2.control(i).W = WC;
end
[wellSolsC2, statesC2] = simulateScheduleAD(state0C1, modelC2, scheduleC2);


%% Compare the well solutions
plotWellSols({wellSols, wellSolsC1, wellSolsC2}, cumsum(schedule.step.val), ...
   'DatasetNames', {'Fine scale', 'Harmonic', 'Flow-based'}, 'Field', 'qOs');

%% Compute flow diagnostics
% As an alternative to looking at well curves, we can also look at the flow
% diagnostics of the models. Flow diagnostics are simple routines based on
% time-of-flight and tracer equations, which aim to give a qualitative
% understanding of reservoir dynamics. Here, we take the state after a
% single time step as a snapshot for both the fine and coarse model and
% compute time-of-flight and well tracers.
mrstModule add diagnostics
D    = computeTOFandTracer(states{2},   G,  rock,  'Wells', schedule.control.W);
DC1  = computeTOFandTracer(statesC1{2}, GC, rockC, 'Wells', scheduleC1.control.W);
DC2  = computeTOFandTracer(statesC2{2}, GC, rockC, 'Wells', scheduleC1.control.W);
WP   = computeWellPairs(states{1},   G,  rock,  schedule.control.W,   D);
WPC1 = computeWellPairs(statesC1{1}, GC, rockC, scheduleC1.control.W, DC1);
WPC2 = computeWellPairs(statesC2{1}, GC, rockC, scheduleC1.control.W, DC2);

figure, plotWellAllocationComparison(DC1, WPC1, D, WP);%, 'plotrow', true);
figure, plotWellAllocationComparison(DC2, WPC2, D, WP);%, 'plotrow', true);

%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.
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
