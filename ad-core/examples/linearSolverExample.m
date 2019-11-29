mrstModule add ad-core ad-blackoil ad-props blackoil-sequential linearsolvers agmg
%% Set up model
if ~exist('n', 'var')
    n = 20;
end
% Enable plotting
if ~exist('doPlot', 'var')
    doPlot = true;
end

if ~exist('solveDirect', 'var')
    solveDirect = n^3 <= 100000;
end
assert(isscalar(n));
assert(n <= 100);
hasAGMG = ~isempty(mrstPath('agmg'));

%%
G = cartGrid([n, n, n], [1000, 1000, 100]);
G = computeGeometry(G);
rock = makeRock(G, 1*darcy, 1);
ncells = G.cells.num;

W = [];
W = verticalWell(W, G, rock, 1, 1, n, 'comp_i', [1, 0, 0], 'val', 100*barsa);
W = verticalWell(W, G, rock, 1, n, 1, 'comp_i', [1, 0, 0], 'val', 50*barsa);
W = verticalWell(W, G, rock, n, 1, 1, 'comp_i', [1, 0, 0], 'val', 150*barsa);
W = verticalWell(W, G, rock, n, n, n, 'comp_i', [1, 0, 0], 'val', 200*barsa);

s0 = rand(G.cells.num, 3);
s0 = bsxfun(@rdivide, s0, sum(s0, 2));
state0 = initResSol(G, 50*barsa, s0);
%%
fluid = initSimpleADIFluid('c', [1e-6, 1e-5, 1e-3]/barsa, 'n', [2, 2, 2], 'rho', [1000, 500, 100]);
model = GenericBlackOilModel(G, rock, fluid);
% model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);

tmodel = TransportModel(model);
pmodel = PressureModel(model);

forces = model.getValidDrivingForces();
forces.W = W;

model = model.validateModel(forces);
pmodel = pmodel.validateModel(forces);
tmodel = tmodel.validateModel(forces);

state0 = model.validateState(state0);
dt = 1*year;
%% Utilities
maxIter = 250;
tol = 1e-6;

getter = @(x, fld) cellfun(@(x) x.(fld), x);
getData = @(x) [getter(x, 'PreparationTime'); getter(x, 'LinearSolutionTime'); getter(x, 'PostProcessTime')]';

%% Get systems
problem = model.getEquations(state0, state0, dt, forces);
problem = problem.assembleSystem();
ncomp = model.getNumberOfComponents();
ordering = getCellMajorReordering(ncells, ncomp);

% Pressure
pproblem = pmodel.getEquations(state0, state0, dt, forces);
% Skip wells
pproblem.equations = pproblem.equations(1);
pproblem.equations{1}.jac = pproblem.equations{1}.jac(1);
pproblem = pproblem.assembleSystem();

% Transport
statep = standaloneSolveAD(state0, pmodel, dt, 'W', W, 'LinearSolver', AMGCLSolverAD('keepNumber', ncells));
statet = tmodel.validateState(statep);
state0t = tmodel.validateState(state0);

tproblem = tmodel.getEquations(state0t, statet, dt, forces);
tproblem = tproblem.assembleSystem();

%% Compare CPR
innerTol = 1e-4;
block_arg = {'variableOrdering', ordering, 'equationOrdering', ordering};
base_arg = {'tolerance', tol, 'maxIterations', maxIter, 'keepNumber', ncomp*ncells};
% base_arg = {}
cpr = CPRSolverAD(base_arg{:});
cpr_agmg = CPRSolverAD(base_arg{:},'ellipticSolver', AGMGSolverAD('tolerance', innerTol));
cpr_amgcl = CPRSolverAD(base_arg{:},'ellipticSolver', AMGCLSolverAD('reuseMode', 2, 'tolerance', innerTol));
cpr_cl = AMGCL_CPRSolverAD(base_arg{:}, block_arg{:});
bl = BackslashSolverAD();
gmilu = GMRES_ILUSolverAD(base_arg{:});
%%
fi_solvers = {};
if solveDirect
    fi_solvers{end+1} = bl;
end
fi_solvers{end+1} = gmilu;
fi_solvers{end+1} = cpr;
if hasAGMG
    fi_solvers{end+1} = cpr_agmg;
end
fi_solvers{end+1} = cpr_amgcl;
fi_solvers{end+1} = cpr_cl;

[descr_fi, names_fi] = cellfun(@getDescription, fi_solvers, 'UniformOutput', false);

ns = numel(fi_solvers);
reports_fi = cell(1, ns);
for sno = 1:ns
    fprintf('Fully-implicit: Solving solver %d of %d: %s\n', sno, ns, names_fi{sno});
    solver = fi_solvers{sno};
    [dx, result, reports_fi{sno}] = solver.solveLinearProblem(problem, model);
end
fprintf('Solve done.\n');
%%
if doPlot
    figure(1); clf;
    plotStuff(reports_fi, names_fi, getData);
    title('Fully-implicit system');
end
%% Compare pressure system
parg = {'tolerance', tol, 'maxIterations', maxIter, 'keepNumber', ncells};
amg_stuben = AMGCLSolverAD(parg{:}, 'coarsening', 'ruge_stuben');
amg_aggr = AMGCLSolverAD(parg{:}, 'coarsening', 'aggregation');
amg_smoothed_aggr = AMGCLSolverAD(parg{:}, 'coarsening', 'smoothed_aggregation');
agmg = AGMGSolverAD(parg{:});

p_solvers = {};
if solveDirect
    p_solvers{end+1} = bl;
end
p_solvers{end+1} = amg_stuben;
p_solvers{end+1} = amg_aggr;
p_solvers{end+1} = amg_smoothed_aggr;
if hasAGMG
    p_solvers{end+1} = agmg;
end

[descr_p, names_p] = cellfun(@getDescription, p_solvers, 'UniformOutput', false);

nsp = numel(p_solvers);
reports_p = cell(1, nsp);

for sno = 1:nsp
    fprintf('Pressure: Solving solver %d of %d\n', sno, nsp);
    solver = p_solvers{sno};
    [dx, result, reports_p{sno}] = solver.solveLinearProblem(pproblem, model);
end
if doPlot
    figure(2); clf;
    plotStuff(reports_p, names_p, getData);
    title('Pressure system');
end
%% Compare transport system
% t-solvers: GMRES ilu, backslash, amgcl ilu0, ilu0-block, iluk-block,
% spai0

topts = ['preconditioner', 'relaxation', block_arg];
amgcl_gs = AMGCLSolverAD(topts{:}, base_arg{:}, 'block_size', 1, 'relaxation', 'gauss_seidel');
amgcl_ilu  = AMGCLSolverAD(topts{:}, base_arg{:}, 'block_size', 1, 'relaxation', 'ilu0');
amgcl_bilu = AMGCLSolverAD(topts{:}, base_arg{:}, 'block_size', ncomp, 'relaxation', 'ilu0');
amgcl_bgs = AMGCLSolverAD(topts{:}, base_arg{:}, 'block_size', ncomp, 'relaxation', 'gauss_seidel');

t_solvers = {};
if solveDirect
    p_solvers{end+1} = bl;
end
t_solvers{end+1} = gmilu;
t_solvers{end+1} = amgcl_gs;
t_solvers{end+1} = amgcl_ilu;
t_solvers{end+1} = amgcl_bgs;
t_solvers{end+1} = amgcl_bilu;
[descr_t, names_t] = cellfun(@getDescription, t_solvers, 'UniformOutput', false);

nst = numel(t_solvers);
reports_t = cell(1, nst);

for sno = 1:nst
    fprintf('Transport: Solving solver %d of %d\n', sno, nst);
    solver = t_solvers{sno};
    [dx, result, reports_t{sno}] = solver.solveLinearProblem(tproblem, model);
end
%%
if doPlot
    figure(3); clf;
    plotStuff(reports_t, names_t, getData);
    title('Transport system')
end
%%
function plotStuff(reports, names, getData)
    timing = getData(reports);
    ok = true(size(timing, 1), 1);
    for i = 1:size(timing, 1)
        if isfield(reports{i}, 'Converged')
            ok(i) = reports{i}.Converged;
        end
    end
    timing = timing.*ok;
    bh = bar(timing, 'stacked');
    legend('Pre', 'Solve', 'Post');
    set(gca, 'XTickLabel', names, 'TickLabelInterpreter', 'none');
    set(gca, 'XTickLabelRotation', -90)
    tot = sum(timing, 2);
    for i = 1:numel(tot)
        r = reports{i};
        if isfield(r, 'Iterations') && r.Iterations(end) > 0
            if ok(i)
                s = num2str(r.Iterations(end));
            else
                s = 'Fail';
            end
            th = text(i, tot(i) + 0.05*max(tot), s);
            set(th, 'HorizontalAlignment', 'center', 'FontSize', 14)
        end
    end
    ylim([0, 1.1*max(tot)]);
end