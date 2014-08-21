mrstModule add ad-fi deckformat ad-core ad-blackoil

% use new code for wells; currently located at
% projects/wells_ad-fi

% mrstVerbose true

% Read and process file.
filedir = '~/simmatlab/projects/ad-fi_benchmarks/b2-SPE9/';
% filedir = 'D:/Jobb/ad-fi_benchmarks/b2-SPE9/';
fn    = fullfile(filedir, 'BENCH_SPE9.DATA');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

fluid = initDeckADIFluid(deck);
gravity on

% Initial solution
p0  = deck.SOLUTION.PRESSURE;
sw0 = deck.SOLUTION.SWAT;
sg0 = deck.SOLUTION.SGAS;
s0  = [sw0, 1-sw0-sg0, sg0];
rs0 = deck.SOLUTION.RS;

state = struct('s', s0, 'rs', rs0, 'pressure', p0);   clear k p0 s0 rs0
state.rv = 0;


clear boModel
clear nonlinear
% clear CPRSolverAD
clear linsolve

model = selectModelFromDeck(G, rock, fluid, deck);
model.drsMaxRel = .2;
model.dpMaxRel  = .2;
model.dsMaxAbs  = .05;

schedule = convertDeckScheduleToMRST(G, model, rock, deck);

%%
% amgsolver = AGMGSolverAD();
amgsolver = BackslashSolverAD();


linsolve = CPRSolverAD('ellipticSolver', amgsolver, 'diagonalTol', .5);
% 
% linsolve = mldivideSolverAD();

% nonlinear = nonlinearSolver();
% [state, status] = nonlinear.solveTimestep(state, 1*day, boModel)

timer = tic();
[wellSols, states] = simulateScheduleAD(state, model, schedule, 'linearSolver', linsolve);
t_class = toc(timer);

%%
compressedSchedule = compressSchedule(schedule);

rampup = 1*day;
timestepper = GustafssonLikeStepSelector('targetIterationCount', 5,...
                                         'minRelativeAdjustment', sqrt(eps),...
                                         'maxRelativeAdjustment', inf, ...
                                         'firstRampupStep',       rampup, ...
                                         'verbose', true);
nonlinear = NonLinearSolver('timeStepSelector', timestepper, 'verbose', true);
timer = tic();
[ws, s, reports] =  simulateScheduleAD(state, model, compressedSchedule, ...
    'nonlinearSolver', nonlinear, 'linearSolver', linsolve, 'OutputMinisteps', true);
t_step = toc(timer);

