%% Two phases flow simulation using AD and including adjoint computations
% The problem is defined in 'INPUT_NUMGRAD.DATA' which is a simple
% $10\times1\times10$ Cartesian grid with uniform permeability. We simulate
% the evolution of the reservoir for a given schedule and, using the
% adjoint equations, we compute the derivative of a given net present value
% (NPV) function with respect to the control variables.

%% Read the problem from a deckfile
% The problem is defined in 'INPUT_NUMGRAD.DATA' which is a simple
% $10\times1\times10$ Cartesian grid with uniform permeability. We read the
% deck and create the grid, rock and fluid structures from the resulting
% output. This requires the deckformat module.

% This requires the deckformat module.
mrstModule add deckformat ad-core ad-fi optimization ad-props

current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'simple10x1x10.data');
deck = readEclipseDeck(fn);

% Convert to MRST units (SI)
deck = convertDeckUnits(deck);

% Create grid
G = initEclipseGrid(deck);

% Set up the rock structure
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

% Create fluid
fluid = initDeckADIFluid(deck);

% Get schedule
schedule = deck.SCHEDULE;

% Enable this to get convergence reports when solving schedules
verbose = false;

%% Visualize the fluid properties (Relative permeability)
% While the geometry is trivial, we can examine the fluid properties. We do
% this by creating a vector containing the whole range of saturation values
% and sampling the functions contained in the fluid object.
s = 0:0.05:1;
clf;
subplot(2,1,1)
plot(s, fluid.krW(s))
title('kr_w')

subplot(2,1,2)
plot(s, fluid.krOW(s))
title('kr_{ow}')

%% Show the schedules
% There are three different well setups in the schedule. We are interested
% in the gradient based on these wells. We convert the timesteps to days to
% get better plots and visualize the controls for all the timesteps to show
% timescales.
inj = vertcat(schedule.control.WCONINJE) %#ok intentional display
prod = vertcat(schedule.control.WCONPROD) %#ok intentional display
timesteps = convertTo(cumsum(schedule.step.val), day);
controls  = schedule.step.control;

% Extract rates and pressures
rates     = convertTo([inj{:, 5}], stb/day);
pressures = convertTo([prod{:, 9}], milli*barsa);

% Plot pressures and rates in each timestep.
clf;
subplot(2,1,1)
plot(timesteps, pressures(controls), '-x')
title('Producer BHP');
ylabel('millibar')
xlabel('days')
axis(axis() + [0 0 -.1*max(pressures) .1*max(pressures)])

subplot(2,1,2)
plot(timesteps, rates(controls), '-x')
title('Injector rate');
ylabel('stb/day')
xlabel('days')

% Adjust the axis a bit
axis(axis() + [0 0 -.1*max(rates) .1*max(rates)])

%% Visualize fluid properties (B)
% B_O and B_W relate surface volumes of the fluids to reservoir conditions.
% Since we are not dealing with a gas phase, the volume ratio between
% surface and reservoir conditions is not very significant.
p = reshape(0.1*barsa : 10*milli*barsa : 1*barsa, [], 1);

clf;
subplot(2,1,1)
plot(p, 1./fluid.bW(p))
title('B_W')
xlabel('Pressure (Pascal)')

subplot(2,1,2)
plot(p, 1./fluid.bO(p))
title('B_O')
xlabel('Pressure (Pascal)')

%% Compute constants
% Once we are happy with the grid and rock setup, we compute
% transmissibilities. For this we first need the centroids.
G = computeGeometry(G);
T = computeTrans(G, rock);

%% Set up reservoir
% We turn on gravity and set up reservoir and scaling factors.
gravity on

state = initResSol(G, deck.PROPS.PVCDO(1), [.15, .85]);

% The scaling factors are use in the Newton solver and hopefully reduce the
% ill-conditionness of the system. The default values are 1.
scalFacs.pressure = 100*barsa;
scalFacs.rate     = 100/day;

%% Run the whole schedule
% We setup the oil water system by calling the function |initADIsystem| and
% the equations are solved implicitely by calling |runScheduleADI|.
timer = tic;
system = initADISystem({'Oil', 'Water'}, G, rock, fluid);
[wellSols, states] = runScheduleADI(state, G, rock, system, schedule);
t_forward = toc(timer);

%% Create objective functions
% We can then create objective functions, which are here net profit worth.
% Since the adjoint formulation uses one forward run (runScheduleADI) to
% get the values for the objective function and one backward run
% (runAdjointADI) to get the gradients, we create an objective function
% based on the earlier solution. Since we will compare with the more
% computationally intensive numerical gradient, we also define an objective
% function which will be used for approximating the gradient of the
% objective function in a difference scheme. We set up approximate prices
% in USD for both the oil price and the injection cost of the different
% phases.

prices = {'OilPrice',            100  , ...
          'WaterProductionCost',   1  , ...
          'WaterInjectionCost',    0.1, ...
          'DiscountFactor',        0.1 };

objective_adjoint = @(tstep)NPVOW(G, wellSols, schedule, 'ComputePartials', true, 'tStep', ...
                                  tstep, prices{:});
objective_numerical = @(wellSols)NPVOW(G, wellSols, schedule, prices{:});

%% Compute derivatives using the adjoint formulation
% We pass in the objective function of the previous run. The objective
% function had a cost equal to one simulation for each timestep. The
% backward simulation requires one simulation per timestep, for a total of
% two simulations per timestep in the schedule to find the gradient. Note
% that the reverse timesteps are much cheaper to compute because only one
% iteration is required per timestep since the system is linear. Because
% the setup in this example is relatively small, we store the states in
% memory and input them as a keyword argument. If the forward simulation
% would be too big to store in full in memory, it can be saved to disk by
% running runScheduleADI with the 'writeOutput' parameter set to true.
% runAdjointADI will then read the files if not supplied with forward
% simulations directly. ('ForwardStates' not supplied or empty)
timer = tic;
adjointGradient = runAdjointADI(G, rock, fluid, schedule, objective_adjoint, system,  'Verbose', verbose, 'ForwardStates', states);
t_adjoint = toc(timer);

%% Find gradients numerically
% To find the numerical gradients we need to solve each timestep three
% times: Once for the baseline value in the actual timestep, and once for
% each of the two wells to compute the derivative of the objective function
% based on that well.
timer = tic;
numericalGradient = computeNumGrad(state, G, rock, system, schedule, objective_numerical, 'scaling', scalFacs, 'Verbose', verbose);
t_gradient = toc(timer);

%% Plot the gradients
% We find the gradient for each well and for each unique well setup in the
% schedule. There are three distinct well setups in the current schedule.
% The gradients are plotted per well with data points corresponding to the
% unique schedules, showing that they are visually indistinguishable.

wellNames = {wellSols{1}.name};

ga = cell2mat(adjointGradient);
gn = cell2mat(numericalGradient);
figure(1)
clf;
subplot(1,2,1)
% The first well is controlled by rate.
plot(ga(1,:)/day,'-ob'), hold on
plot(gn(1,:)/day,'-xr')
title(['Well 1 (', wellNames{1}, ')'])
xlabel('Control #')

subplot(1,2,2)
% The second well is controlled by pressure.
plot(ga(2,:)*barsa,'-ob'), hold on
plot(gn(2,:)*barsa,'-xr'), hold on
title(['Well 2 (', wellNames{2}, ')'])
xlabel('Control #')
legend({'Adjoint', 'Numerical'})

%% Plot the time taken
% Since the adjoint formulation uses both a forward and a backward
% simulation, these must be plotted together.
figure(2)
clf;
bar([t_forward, t_adjoint, 0; 0, 0, t_gradient], 'barlayout', 'stacked');
set(gca, 'XTickLabel', {'Adjoint', 'Numerical'})
legend('Adjoint (Forward)', 'Adjoint (Backwards)', 'Numerical', 'location', 'NorthOutside')
ylabel('t (seconds)')
title('Time used')

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
