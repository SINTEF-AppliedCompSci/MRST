%% Example: Depletion of a closed or open reservoir compartment
% In this tutorial, we will show how to set up a simulator from scratch in
% the automatic differentiation, object-oriented (AD-OO) framework without
% the use of input files. As an example we consider a 2D rectangular
% reservoir compartment with homogeneous properties, drained by a single
% producer at the midpoint of the top edge. The compartment is either
% closed (i.e., sealed by no-flow boundary conditions along all edges), or
% open with constant pressure support from an underlying, infinite aquifer,
% which we model as a constant-pressure boundary condition.

mrstModule add ad-props  ad-core ad-blackoil

%% Grid, petrophysics, and fluid objects
% To create a complete model object using the AD-OO framework, we first
% need to define three standard MRST structures representing the grid and
% the rock and fluid properties

% The grid and rock model
G    = computeGeometry(cartGrid([50 50],[1000 100]));
rock = makeRock(G, 100*milli*darcy, 0.3);

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
gravity reset on
wModel = WaterModel(G, rock, fluid,'gravity',[0 -norm(gravity)]);
% Prepare the model for simulation.
wModel = wModel.validateModel();
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
wc = sub2ind(G.cartDims, floor(G.cartDims(1)/2), G.cartDims(2));
W = addWell([], G, rock,  wc,     ...
        'Type', 'bhp', 'Val', pR-50*barsa, ...
        'Radius', 0.1, 'Name', 'P1','Comp_i',1,'sign',1);

% Boundary conditions: fixed pressure at bottom and no-flow elsewhere
bc=pside([],G,'South',200*barsa,'sat',1);

% Schedule: describing time intervals and corresponding drive mechanisms
schedule1 = simpleSchedule(diff(linspace(0,5*day,41)), 'bc', bc, 'W', W);
schedule2 = simpleSchedule(diff(linspace(0,5*day,41)), 'W', W);

%% Reservoir state
% The last component we need in order to specify our reservoir model is the
% reservoir state, i.e., the fluid pressure. For multiphase models, the
% state also includes the phase saturations and compositions. In our case,
% we first set a constant pressure, and call on a solver from |ad-core| to
% compute vertical equilibrium.
state = initResSol(G, pR); % Constant pressure
state = wModel.validateState(state);

% Vertical equilibrium
verbose = false;
nonlinear = NonLinearSolver();
state = nonlinear.solveTimestep(state, 10000*day, wModel, 'bc', bc);

clf,
plotCellData(G,state.pressure/barsa,'EdgeColor','none');
colorbar

%% Run simulations
% To make the case a bit more interesting, we compute the flow problem
% twice. The first simulation uses the prescribed boundary conditions,
% which will enable fluids to pass out of the north boundary. In the second
% simulation, we close the system by imposing no-flow conditions also on
% the north boundary

% Simulation pressure 200 bar at top
[wellSols1, states1] = simulateScheduleAD(state, wModel, schedule1);

% Simulation with no-flow at top
[wellSols2, states2] = simulateScheduleAD(state, wModel, schedule2);

%% Plot results
% Prepare animation
wpos = G.cells.centroids(wc,:); clf
set(gcf,'Position',[600 400 800 400])

h1 = subplot('Position',[.1 .11 .34 .815]);
hp1 = plotCellData(G,state.pressure/barsa,'EdgeColor','none');
title('Open compartment w/aquifer'); caxis([pR/barsa-50 pR/barsa]);

h2 = subplot('Position',[.54 .11 .4213 .815]);
hp2 = plotCellData(G,state.pressure/barsa,'EdgeColor','none');
title('Closed compartment'); caxis([pR/barsa-50 pR/barsa]);

colormap(jet(32)); colorbar

% Animate solutions by resetting CData property of graphics handle
for i=1:numel(states1)
    set(hp1,'CData', states1{i}.pressure/barsa);
    set(hp2,'CData', states2{i}.pressure/barsa);
    drawnow, pause(0.1);
end

% Launch plotting of well responses
plotWellSols({wellSols1,wellSols2}, ...
    'Datasetnames',{'Aquifer','Closed'}, 'field','qWr');

%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
