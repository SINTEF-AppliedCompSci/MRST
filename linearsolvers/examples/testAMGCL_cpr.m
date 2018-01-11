%% Example demonstrating AMGCL on a few test problems
mrstModule add msrsb incomp coarsegrid spe10 linearsolvers ad-core ad-blackoil ad-props
%% Setup problem

nx = 1000;
ny = nx;
G = cartGrid([nx, ny, 1], [1, 1, 1]);
G = computeGeometry(G);

rock = makeRock(G, 0.1*darcy, 0.5);

fluid = initSimpleADIFluid('c', [1, 1, 1]*1e-5/barsa);
model = TwoPhaseOilWaterModel(G, rock, fluid);

state0 = initResSol(G, 100*barsa, [0.5, 0.5]);

bc = [];
src = [];
W = [];

bc = pside(bc, G, 'xmin', 50*barsa, 'sat', [1, 0]);
bc = pside(bc, G, 'xmax', 50*barsa, 'sat', [1, 0]);

forces = struct('bc', bc, 'W', W, 'src', src);
model = model.validateModel(forces);
state0 = model.validateState(state0);


dt = 30*day;
state = state0;
[problem, state] = model.getEquations(state0, state, dt, forces);

problem = problem.assembleSystem();


A0 = problem.A;
b0 = problem.b;

lsolve = BackslashSolverAD();
lsolve.keepNumber = (model.water + model.gas + model.oil)*G.cells.num;
[A0, b0, sys] = lsolve.reduceLinearSystem(A0, b0);

ncomp = model.water + model.oil + model.gas;

subs = (1:ncomp*G.cells.num);
subs = reshape(subs, [], ncomp)';
subs = subs(:);
b0 = b0./norm(b0, inf);

A = A0(subs, subs);
b = b0(subs);

%% 
clc
At = A';
tic();
[x0, err] = amgcl_matlab_cpr(At, b, 1e-6, 10, 1, 1, 1, 1, 2);
t_amg = toc();
toc()

% ref = A\b;

err

%%
% ref = A\b;

tic();
[x1, err] = callAMGCL_cpr(At, b, 2, 'isTransposed', true, 'maxIterations', 10, 'cellMajorOrder', true);
t_wrapper = toc();
% norm(ref - x1)/norm(ref)

%%
% ref = A0\b0;

tic();
[x2, err] = callAMGCL_cpr(A0, b0, 2, 'isTransposed', false, 'maxIterations', 10, 'cellMajorOrder', false);
t_inline = toc();
% norm(ref - x2)/norm(ref)

%%
% !rm /home/moyner/bitbucket/mrst-solvers/linearsolvers/amgcl/utils/amgcl_matlab_cpr.mexa64

%%
[G, rock, fluid, deck, state] = setupSPE9();
model = selectModelFromDeck(G, rock, fluid, deck);
schedule = convertDeckScheduleToMRST(model, deck);
%%
lsolve = AMGCL_CPRSolverAD('block_size', 3, 'maxIterations', 100, 'tolerance', 1e-6);
% lsolve.applyRightDiagonalScaling = true;
% lsolve.solver = 'gmres';
% lsolve.t_relaxation = 'ilu0';

[wellSols, states, report] = simulateScheduleAD(state, model, schedule, 'linearsolver', lsolve);

%%
