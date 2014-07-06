%% Read the problem from a deckfile
% The problem is defined in 'INPUT_NUMGRAD.DATA' which is a simple
% $10\times1\times10$ Cartesian grid with uniform permeability. We read the
% deck and create the grid, rock and fluid structures from the resulting
% output. This requires the deckformat module.
mrstModule add deckformat ad-fi ad-refactor

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

% Oil rel-perm from 2p OW system.
% Needed by equation implementation function 'eqsfiOWExplictWells'.
fluid.krO = fluid.krOW;

% Get schedule
schedule = deck.SCHEDULE;

% Enable this to get convergence reports when solving schedules
verbose = false;

%% Compute constants
% Once we are happy with the grid and rock setup, we compute
% transmissibilities. For this we first need the centroids.
G = computeGeometry(G);
T = computeTrans(G, rock);

%% Set up reservoir
% We turn on gravity and set up reservoir and scaling factors.
gravity on

state = initResSol(G, deck.PROPS.PVCDO(1), [.15, .85]);

scalFacs.pressure = 100*barsa;
scalFacs.rate     = 100/day;
%% Run the whole schedule
% This is done to get values for the wells for all timesteps. Since the
% case is fairly small,
timer = tic;
system = initADISystem({'Oil', 'Water'}, G, rock, fluid, 'cpr', true);
system.nonlinear.cprBlockInvert = false;


[wellSols states] = runScheduleADI(state, G, rock, system, schedule);
t_forward = toc(timer);


%%
for i = 1:numel(schedule.control)
    W = processWellsLocal(G, rock, schedule.control(i), ...
                                     'Verbose', false, ...
                                     'DepthReorder', false);
    schedule.control(i).W = W;
end

%% New way...
clear owModel
clear nonlinear

% multiscaleSolver = multiscaleVolumeSolverAD(CG);

owModel = TwoPhaseOilWaterModel(G, rock, fluid, 'deck', deck);
% linsolve = CPRSolverAD('ellipticSolver', multiscaleSolver);
linsolve = CPRSolverAD();
% linsolve = mldivideSolverAD();
% nonlinear = nonlinearSolver();
% [state, status] = nonlinear.solveTimestep(state, 1*day, boModel)


[wellSols, states] = simulateScheduleAD(state, owModel, schedule, 'linearSolver', linsolve);