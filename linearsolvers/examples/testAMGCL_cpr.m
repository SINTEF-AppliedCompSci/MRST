%% Example demonstrating AMGCL on a few test problems
mrstModule add msrsb incomp coarsegrid spe10 linearsolvers ad-core ad-blackoil ad-props
%% Setup problem

nx = 100;
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


A = problem.A;
b = problem.b;

lsolve = BackslashSolverAD();
lsolve.keepNumber = (model.water + model.gas + model.oil)*G.cells.num;
[A, b, sys] = lsolve.reduceLinearSystem(A, b);

ncomp = model.water + model.oil + model.gas;

subs = (1:ncomp*G.cells.num);
subs = reshape(subs, [], ncomp)';
subs = subs(:);


A = A(subs, subs);
b = b(subs);

%% 
At = A';
tic();
[x, err] = amgcl_matlab_cpr(At, b, 1e-6, 10, 1, 1, 1, 1, 2);
t_amg = toc();

% ref = A\b;

err

%%
!rm /home/moyner/bitbucket/mrst-solvers/linearsolvers/amgcl/utils/amgcl_matlab_cpr.mexa64