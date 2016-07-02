%% Example: Single-Phase Water Flow
% In this tutorial, we will show how to set up a simulator from scratch in
% the automatic differentiation, object-oriented (AD-OO) framework without
% the use of input files. As an example we consider a 2D rectangular
% reservoir with homogeneous properties, with a single well at the midpoint
% of the south edge and fixed pressure prescribed along the north edge.

mrstModule add ad-props  ad-core ad-blackoil

%% Grid, petrophysics, and fluid objects
% To create a complete model object using the AD-OO framework, we first
% need to define three standard MRST structures representing the grid and
% the rock and fluid properties

% The grid
G = cartGrid([50 50],[1000 100]);
G = computeGeometry(G);

% Permeability and porosity
rock.perm  = 100*milli*ones(G.cells.num,1)*darcy;
rock.poro  = ones(G.cells.num,1)*0.3;

% Transmissibility
T  = computeTrans(G, rock);

% Fluid properties
pR  = 200*barsa;
fluid = initSimpleADIFluid('phases','W',           ... % Fluid phase: water
                           'mu',  1*centi*poise,   ... % Viscosity
                           'rho', 1000,            ... % Surface density [kg/m^3]
                           'c',   1e-4/barsa,      ... % Fluid compressibility
                           'cR',  1e-5/barsa       ... % Rock compressibility
                           );

%% Make Reservoir Model
% We can now use the three objects defined above to instantiate the
% reservoir model. To this end, we will use a WaterModel, which is a
% specialization of the general ReservoirModel implemented in |ad-core|.
% The only extra thing we need to do is to explicitly set the gravity
% direction. By default, the gravity in MRST is a 3-component vector that
% points in the positive z-direction. Here, we set it to a 2-component
% vector pointing in the negative y-direction.
gravity('y'); grav = gravity;
wModel = WaterModel(G, rock, fluid,'gravity',-grav(1:2));

% Enable this to get convergence reports when solving schedules

%% Drive mechansims and schedule
% The second thing we need to specify is the mechanisms that will drive
% flow in the reservoir, i.e., the wells and boundary conditions. These may
% change over time and MRST therefore uses the concept of a schedule that
% describes how the drive mechansims change over time. In our case, we use
% the same setup for the whole simulation. The schedule also enables us to
% specify the time steps we wish the simulator to use, provided that these
% give convergent and stable computations. (If not, the simulator may cut
% the time step).

% Well: at the midpoint of the south edge
wc = sub2ind( G.cartDims, floor(G.cartDims(1)/2), 1);
W = addWell([], G, rock,  wc,     ...
        'Type', 'bhp', 'Val', 100*barsa+pR, ...
        'Radius', 0.1, 'Name', 'P1','Comp_i',1,'sign',1);

% Boundary conditions: fixed pressure at top and no-flow elsewhere
bc=pside([],G,'North',200*barsa,'sat',1);

% Schedule: describing time intervals and corresponding drive mechanisms
dt = diff(linspace(0,1*day,20));
schedule = struct(...
    'control', struct('W',{W}), ...
    'step',    struct('control', ones(numel(dt),1), 'val', dt));
for i=1:numel(schedule.control)
    schedule.control(i).bc = bc;
    for j=1:numel(schedule.control.W)
        schedule.control(i).W(j).compi=1;        % Saturation in wellbore
    end
end

%% Reservoir state
% The last component we need in order to specify our reservoir model is the
% reservoir state, i.e., the fluid pressure. For multiphase models, the
% state also includes the phase saturations and compositions. In our case,
% we first set a constant pressure, and call on a solver form |ad-core| to
% compute vertical equilibrium.
state.pressure = ones(G.cells.num,1)*pR;                % Constant pressure
state.wellSol  = initWellSolAD([], wModel, state);      % No well initially

% Vertical equilibrium
verbose = false;
nonlinear = NonLinearSolver;
state = nonlinear.solveTimestep(state, 10000*day, wModel, 'bc', bc);

clf,
plotCellData(G,state.pressure/barsa,'EdgeColor','none');
colorbar

%% Run Simulations and Plot Results
% To make the case a bit more interesting, we compute the flow problem
% twice. The first simulation uses the prescribed boundary conditions,
% which will enable fluids to pass out of the north boundary. In the second
% simulation, we close the system by imposing no-flow conditions also on
% the north boundary

% Simulation pressure 200 bar at top
[wellSols1, states1] = simulateScheduleAD(state, wModel, schedule);

% Simulation with no-flow at top
schedule2 = schedule;
for i=1:numel(schedule2.control), schedule2.control(i).bc = []; end
[wellSols2, states2] = simulateScheduleAD(state, wModel, schedule2);

% Prepare animation
wpos = G.cells.centroids(wc,:); clf
set(gcf,'Position',[600 400 800 400]); colormap(jet(32));
h1 = subplot('Position',[.1 .11 .34 .815]);
title('With pressure b.c'); caxis([200 300]); hp1 = [];
h2 = subplot('Position',[.54 .11 .4213 .815]);
title('With no-flow b.c.'); caxis([200 300]); hp2 = [];
colorbar

% Animate solutions
for i=1:numel(states1)
    subplot(h1); delete(hp1);
    hp1 = plotCellData(G,states1{i}.pressure/barsa,'EdgeColor','none');
    %hold on, plot(wpos(1),wpos(2),'o','Color','r','MarkerSize',10), hold off
    % caxis([200 300])
    
    subplot(h2); delete(hp2);
    plotCellData(G,states2{i}.pressure/barsa,'EdgeColor','none');
    %hold on, plot(wpos(1),wpos(2),'o','Color','r','MarkerSize',10), hold off
    % caxis([200 300])
    % colorbar
    
    drawnow;
end

% Launch plotting of well responses
plotWellSols({wellSols1,wellSols2}, ...
    'Datasetnames',{'Pressure','No-flow'}, 'field','qWr');

%%
