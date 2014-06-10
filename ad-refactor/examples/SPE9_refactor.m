mrstModule add ad-fi deckformat

% use new code for wells; currently located at
% projects/wells_ad-fi

mrstVerbose true

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

schedule = deck.SCHEDULE;
system = initADISystem(deck, G, rock, fluid, 'cpr', true);
% Hack until this is handled by initADISystem
system.activeComponents.vapoil = 0;
% switch off individual well solves:
system.stepOptions.solveWellEqs = false;
% use new cpr based on dnr
system.nonlinear.cprBlockInvert = false;
% convergence is overall better for quite strict limits on update 
system.stepOptions.drsMax = .2;
system.stepOptions.dpMax  = .2;
system.stepOptions.dsMax  = .05;
% gmres tol needs to be quite strict
system.nonlinear.cprRelTol = 1e-3;
system.pscale = 1/(200*barsa);

schedule.step.control = schedule.step.control(1:10);
schedule.step.val = schedule.step.val(1:10);

%%
timer = tic;
[wellSols, states, iter] = runScheduleADI(state, G, rock, system, schedule);
t_standard = toc(timer);


%%
for i = 1:numel(schedule.control)
    W = processWellsLocal(G, rock, schedule.control(i), ...
                                     'Verbose', false, ...
                                     'DepthReorder', false);
    schedule.control(i).W = W;
end
clear boModel
clear nonlinear
% clear CPRSolverAD
clear linsolve

boModel = threePhaseBlackOilModel(G, rock, fluid, ...
                                        'drsMax', .2,...
                                        'dpMax', .2', ...
                                        'dsMax', .05, ...
                                        'disgas', true, ...
                                        'deck', deck);
                                    
%%
amgsolver = AGMGSolverAD();
% amgsolver = mldivideSolverAD();


linsolve = CPRSolverAD('ellipticSolver', amgsolver);
% 
% linsolve = mldivideSolverAD();

% nonlinear = nonlinearSolver();
% [state, status] = nonlinear.solveTimestep(state, 1*day, boModel)

timer = tic();
[wellSols, states] = runScheduleRefactor(state, boModel, schedule, 'linearSolver', linsolve);
t_class = toc(timer);

%%
cdims = ceil(G.cartDims./[10 10 5]);
p = partitionUI(G, cdims);
% p = partitionMETIS(G, computeTrans(G, rock), prod(cdims));
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);


multiscaleSolver = multiscaleVolumeSolverAD(CG);
linsolvems = CPRSolverAD('ellipticSolver', multiscaleSolver);


timer = tic();
[wellSols, states] = runScheduleRefactor(state, boModel, schedule, 'linearSolver', linsolvems);
t_ms = toc(timer);

