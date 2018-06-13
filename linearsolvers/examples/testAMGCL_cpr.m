%% Example demonstrating AMGCL on a few test problems
mrstModule add msrsb incomp coarsegrid spe10 linearsolvers ad-core ad-blackoil ad-props
%% Setup problem
if 1
    nx = 10;
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
    dt = 30*day;

elseif 1
    mrstModule add spe10
    [state0, model, schedule]  = setupSPE10_AD('layers', 1:5);
    forces = schedule.control(1);
    dt = schedule.step.val(1);
    G = model.G;
else
    [G, rock, fluid, deck, state0] = setupSPE9();
    model = selectModelFromDeck(G, rock, fluid, deck);
    schedule = convertDeckScheduleToMRST(model, deck);
    forces = schedule.control(1);
    dt = schedule.step.val(1);
end
model = model.validateModel(forces);
state0 = model.validateState(state0);


state = state0;
[problem, state] = model.getEquations(state0, state, dt, forces);

problem = problem.assembleSystem();


A0 = problem.A;
b0 = problem.b;

lsolve = BackslashSolverAD();
lsolve.keepNumber = (model.water + model.gas + model.oil)*G.cells.num;
[A0, b0, sys] = lsolve.reduceLinearSystem(A0, b0);


ncomp = model.water + model.oil + model.gas;

pix = 1:G.cells.num;
for i = 2:ncomp
    ix = (1:G.cells.num) + (i-1)*G.cells.num;
    A0(pix, :) = A0(pix, :) + A0(ix, :);
    b0(pix) = b0(pix) + b0(ix);
end

subs = getCellMajorReordering(G.cells.num, ncomp);

b0 = b0./norm(b0, inf);

A = A0(subs, subs);
b = b0(subs);
% b = b0;
its = 100;
At = A';

%%
% ref = A\b;
tic();
[x1, err1] = callAMGCL_cpr(At, b, ncomp, 'isTransposed', true, 'maxIterations', its, 'cellMajorOrder', true);
t_wrapper = toc();  
err1
% norm(ref - x1)/norm(ref)

%%
% ref = A0\b0;

tic();
[x2, err2] = callAMGCL_cpr(A0, b0, ncomp, 'isTransposed', false, 'maxIterations', its, 'cellMajorOrder', false);
t_inline = toc();
err2
% norm(ref - x2)/norm(ref)

%%
% !rm /home/moyner/bitbucket/mrst-solvers/linearsolvers/amgcl/utils/amgcl_matlab_cpr.mexa64

%%
[G, rock, fluid, deck, state] = setupSPE9();
model = selectModelFromDeck(G, rock, fluid, deck);
model.AutoDiffBackend = DiagonalAutoDiffBackend();
model.FacilityModel = UniformFacilityModel(model);
schedule = convertDeckScheduleToMRST(model, deck);
%%
ncomp = model.water+model.gas+model.oil;

lsolve = AMGCLSolverAD('maxIterations', 200, 'tolerance', 1e-3,...
    'preconditioner', 'relaxation',...
    'solver', 'idrs', ...
    'relaxation', 'ilu0');
if 1
    ord = getGridSYMRCMOrdering(model.G);
else
    ord = [];
end

ndof = ncomp*G.cells.num;
ordering = getCellMajorReordering(G.cells.num, ncomp, 'ndof', ndof);
lsolve.variableOrdering = ordering;
lsolve.equationOrdering = ordering;

if 0
    lsolve = AMGCL_CPRSolverAD('block_size', ncomp, 'maxIterations', 200, 'tolerance', 1e-3);
    lsolve.setSRelaxation('ilu0');
    lsolve.trueIMPES = true;
    lsolve.setCoarsening('aggregation');
end
% lsolve.applyRightDiagonalScaling = true;
lsolve.setSolver('bicgstab');

if isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend')
    lsolve.reduceToCell = false;
    lsolve.keepNumber = G.cells.num*ncomp;
end

[wellSols, states, report] = simulateScheduleAD(state, model, schedule, 'linearsolver', lsolve);

%%
ls_time = 0;
for i = 1:numel(report.ControlstepReports)
    rr = report.ControlstepReports{i};
    for j = 1:numel(rr.StepReports{1}.NonlinearReport)
        rrr = rr.StepReports{1}.NonlinearReport{j};
        if rrr.Converged
            continue
        end
        ls_time = ls_time + rrr.LinearSolver.SolverTime;
    end
end
