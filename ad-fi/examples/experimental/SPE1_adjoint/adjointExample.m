%% Read case from file
% SPE1 model

mrstModule add deckformat ad-fi optimization

% Read and process file.
current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'odeh_adi_simplified.data');
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
gravity on


% The initial state is provided as a binary file. The initial state
% contains a uniform mixture of water (.12) and oil (.88).
load initialState;

%% Initialize schedule and system before solving for all timesteps
% We extract the schedule from the read deck and create a ADI system for
% our problem. The system autodetects a black oil problem and sets up
% default values for the various options. The only thing we change is to
% disable the CPR preconditioner as the problem is too small to benefit
% from preconditioning: The overhead required for the preconditioner is
% bigger than the benefits offered by a specialized solver.
%
% During some time steps (67 and 91) the Newton iterations oscillate. The
% solver detects this, and dampens or relaxes the step length when this
% behavior is observed.
%
% To see detailed convergence analysis during each time step, set verbose
% to on using

mrstVerbose on

schedule = deck.SCHEDULE;
system   = initADISystem(deck, G, rock, fluid);

system.nonlinear.maxIterations = 30;

% no relaxation
system.nonlinear.relaxMax    = 0.5;   % default = 0.5
system.nonlinear.relaxRelTol = 0.1;   % default = .2
system.nonlinear.cpr         = true;  % default = false
system.nonlinear.cprRelTol   = 1e-6;  % default = 1e-3
system.nonlinear.linesearch  = false; % default = false
system.nonlinear.lineIter    = 10;    % default = 10

% system.nonlinear.cprEllipticSolver = @(A,b) agmg(A, b, [], [], 100, 0);
system.nonlinear.cprEllipticSolver = @mldivide; % default = []
system.nonlinear.bhpcontrols       = false;     % default = false
system.nonlinear.changeWells       = false;     % default = false
system.nonlinear.cprBlockInvert    = false;     % default = true

% Newton step options

system.stepOptions.dsMax = 0.2; % default = 0.2
% system.stepOptions.dpMax = 1e7; % default = inf
% system.stepOptions.drsMax = inf; % default = inf
system.stepOptions.solveWellEqs = true; % default = true

outputDir = fullfile(ROOTDIR, 'modules', 'ad-fi', 'wells', 'cache2');

system.stepFunction = @(state0, state, meta, dt, W, G, system) ...
   stepBlackOil        (state0, state, meta, dt, G, W, system, fluid);

timer = tic;
[wellSols, states, iter] = ...
   runScheduleADI(state, G, rock, system, schedule, ...
                  'writeOutput', true, 'outputDir', outputDir);
toc(timer)

%% Objective functions
% Create objective functions for the different systems. We set up
% approximate prices in USD for both the oil price and the injection cost
% of the different phases.

prices = { 'OilPrice',            100 , ...
           'GasPrice',            1.  , ...
           'GasInjectionCost',    0.01, ...
           'WaterProductionCost', 1   , ...
           'WaterInjectionCost',  0.1 , ...
           'DiscountFactor',      0.1 };

% We define
objectiveBO = @(tstep) ...
   NPVBlackOil(G, wellSols, schedule, 'ComputePartials', true, ...
               'tStep', tstep, prices{:});

objective = ...
   NPVBlackOil(G, wellSols, schedule, 'ComputePartials', false, prices{:});

%% Compute gradient using the adjoint formulation
% We pass a function handle to the polymer equations and calculate the
% gradient with regards to our control variables. The control variables are
% defined as the last two variables, i.e. well closure (rate/BHP) and
% polymer injection rate.

% Gradient is computed along all wells control (qWs, qOs, qGs, pBHP)
ctrl = [];
adjointGradient = ...
   runAdjointADI(G, rock, fluid, schedule, objectiveBO, system, ...
                 'Verbose', true, 'ForwardStates', states);

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
