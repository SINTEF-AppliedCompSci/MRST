mrstModule add ad-core ad-blackoil ad-props sequential linearsolvers agmg
%% Set up model
if ~exist('n', 'var')
    n = 40;
end
% Enable plotting
if ~exist('doPlot', 'var')
    doPlot = true;
end

if ~exist('solveDirect', 'var')
    solveDirect = true;
end

solveDirect = solveDirect && n^3 <= 50000;
assert(isscalar(n));
assert(n <= 100);
hasAGMG = ~isempty(mrstPath('agmg'));
n = ceil(n);
% Plot offset
offset = 0;
%% Set up test problem
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
%% Set up model
fluid = initSimpleADIFluid('c', [1e-6, 1e-5, 1e-3]/barsa, 'n', [2, 2, 2], 'rho', [1000, 500, 100]);
model = GenericBlackOilModel(G, rock, fluid);
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, ...
                                                'deferredAssembly', true, ...
                                                'rowMajor', true);

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
maxIter = 100;
tol = 1e-6;


%% Get systems
problem = model.getEquations(state0, state0, dt, forces);
ncomp = model.getNumberOfComponents();
ordering = getCellMajorReordering(ncells, ncomp);

% Pressure
pproblem = pmodel.getEquations(state0, state0, dt, forces);
% Skip wells
% pproblem.equations = pproblem.equations(1);
% pproblem.equations{1}.jac = pproblem.equations{1}.jac(1);
% We can pre-assemble to save a bit of time
pproblem = pproblem.assembleSystem();

% Transport
statep = standaloneSolveAD(state0, pmodel, dt, 'W', W, 'LinearSolver', AMGCLSolverAD('keepNumber', ncells));
statet = tmodel.validateState(statep);
state0t = tmodel.validateState(state0);

tproblem = tmodel.getEquations(state0t, statet, dt, forces);
tproblem = tproblem.assembleSystem();

%% Compare fully-implicit
innerTol = 1e-4;
block_arg = {'variableOrdering', ordering, 'equationOrdering', ordering};
base_arg = {'tolerance', tol, 'maxIterations', maxIter, 'keepNumber', ncomp*ncells};
cpr_matlab_arg = [base_arg, 'diagonalTol', nan];
cpr = CPRSolverAD(cpr_matlab_arg{:});
if hasAGMG
    cpr_agmg = CPRSolverAD(cpr_matlab_arg{:},'ellipticSolver', AGMGSolverAD('tolerance', innerTol));
end
cpr_amgcl = CPRSolverAD(cpr_matlab_arg{:},'ellipticSolver', AMGCLSolverAD('reuseMode', 2, 'tolerance', innerTol));
cpr_cl = AMGCL_CPRSolverAD(base_arg{:}, block_arg{:});
cpr_cl_block_w = AMGCL_CPRSolverBlockAD(base_arg{:}, ...
    'coarsening', 'smoothed_aggregation', 'relaxation', 'spai0', ...
    'aggr_eps_strong', 0.1, 'aggr_over_interp', 1.5, ...
    'solver', 'bicgstab', 'npre', 1, 'npost', 2, 'ncycle', 1, 'id', '-bcsr-tweaked');
cpr_cl_block = AMGCL_CPRSolverBlockAD(base_arg{:}, 'id', '-bcsr');

%%
bl = BackslashSolverAD();
gmilu = GMRES_ILUSolverAD(base_arg{:});
% Store
fi_solvers = {};
if solveDirect
    fi_solvers{end+1} = bl;
end
fi_solvers{end+1} = gmilu;
if solveDirect
    fi_solvers{end+1} = cpr;
end
if hasAGMG
    fi_solvers{end+1} = cpr_agmg;
end
fi_solvers{end+1} = cpr_amgcl;
fi_solvers{end+1} = cpr_cl;
fi_solvers{end+1} = cpr_cl_block;
fi_solvers{end+1} = cpr_cl_block_w;

[descr_fi, names_fi] = cellfun(@getDescription, fi_solvers, 'UniformOutput', false);
% Solving
ns = numel(fi_solvers);
reports_fi = cell(1, ns);
for sno = 1:ns
    fprintf('Fully-implicit: Solving solver %d of %d: %s\n', sno, ns, names_fi{sno});
    solver = fi_solvers{sno};
    [dx, result, reports_fi{sno}] = solver.solveLinearProblem(problem, model);
end
fprintf('Solve done.\n');
%% Plot
[timing_fi, its_fi] = getTiming(reports_fi);
if doPlot
    figure(1 + offset); clf;
    plotLinearTimingsForExample(names_fi, timing_fi, its_fi);
    title('Fully-implicit system');
end
%% Setup pressure solver
parg = {'tolerance', tol, 'maxIterations', maxIter, 'keepNumber', ncells};
amg_stuben = AMGCLSolverAD(parg{:}, 'coarsening', 'ruge_stuben', 'id', '-classical');
amg_aggr = AMGCLSolverAD(parg{:}, 'coarsening', 'aggregation', 'id', '-aggregation');
amg_smoothed_aggr = AMGCLSolverAD(parg{:}, 'coarsening', 'smoothed_aggregation', 'id', '-smoothed_aggregation');
amg_smoothed_aggre = AMGCLSolverAD(parg{:}, 'coarsening', 'smoothed_aggr_emin', 'id', '-smoothed_aggr_emin');
if hasAGMG
    agmg = AGMGSolverAD(parg{:});
end

p_solvers = {};
if solveDirect
    p_solvers{end+1} = bl;
end
p_solvers{end+1} = amg_stuben;
p_solvers{end+1} = amg_aggr;
p_solvers{end+1} = amg_smoothed_aggr;
p_solvers{end+1} = amg_smoothed_aggre;

if hasAGMG
    p_solvers{end+1} = agmg;
end

[descr_p, names_p] = cellfun(@getDescription, p_solvers, 'UniformOutput', false);
nsp = numel(p_solvers);
reports_p = cell(1, nsp);

for sno = 1:nsp
    fprintf('Pressure: Solving solver %d of %d: %s\n', sno, nsp, names_p{sno});
    solver = p_solvers{sno};
    [dx, result, reports_p{sno}] = solver.solveLinearProblem(pproblem, model);
end
%%
[timing_p, its_p] = getTiming(reports_p);
if doPlot
    figure(2 + offset); clf;
    plotLinearTimingsForExample(names_p, timing_p, its_p);
    title('Pressure system');
end

%% Compare transport system
% t-solvers: GMRES ilu, backslash, amgcl ilu0, ilu0-block, iluk-block,
% spai0

topts = ['preconditioner', 'relaxation', block_arg];
amgcl_gs = AMGCLSolverAD(topts{:}, base_arg{:}, 'block_size', 1, 'relaxation', 'gauss_seidel', 'id', '-gs');
amgcl_ilu  = AMGCLSolverAD(topts{:}, base_arg{:}, 'block_size', 1, 'relaxation', 'ilu0', 'id', '-ilu0');
amgcl_bilu = AMGCLSolverAD(topts{:}, base_arg{:}, 'block_size', ncomp, 'relaxation', 'ilu0', 'id', '-ilu0');
amgcl_bgs = AMGCLSolverAD(topts{:}, base_arg{:}, 'block_size', ncomp, 'relaxation', 'gauss_seidel', 'id', '-gs');
bamgcl_ilu = AMGCLSolverBlockAD(topts{:}, base_arg{:}, 'block_size', ncomp, 'relaxation', 'ilu0', 'id', '-bcsr-ilu0');

t_solvers = {};
if solveDirect
    t_solvers{end+1} = bl;
end
t_solvers{end+1} = gmilu;
t_solvers{end+1} = amgcl_gs;
t_solvers{end+1} = amgcl_ilu;
t_solvers{end+1} = amgcl_bgs;
t_solvers{end+1} = amgcl_bilu;
t_solvers{end+1} = bamgcl_ilu;

[descr_t, names_t] = cellfun(@getDescription, t_solvers, 'UniformOutput', false);

nst = numel(t_solvers);
reports_t = cell(1, nst);

for sno = 1:nst
    fprintf('Transport: Solving solver %d of %d: %s\n', sno, nst, names_t{sno});
    solver = t_solvers{sno};
    [dx, result, reports_t{sno}] = solver.solveLinearProblem(tproblem, model);
end
%%
[timing_t, its_t] = getTiming(reports_t);

if doPlot
    figure(3 + offset); clf;
    plotLinearTimingsForExample(names_t, timing_t, its_t);
    title('Transport system')
end
%%
function [timing, its] = getTiming(reports)
    timing = getData(reports);
    ok = true(size(timing, 1), 1);
    for i = 1:size(timing, 1)
        if isfield(reports{i}, 'Converged')
            ok(i) = reports{i}.Converged;
        end
    end
    timing = timing.*ok;
    its = zeros(size(reports));
    for i = 1:numel(reports)
        r = reports{i};
        if isfield(r, 'Iterations')
            if ok(i)
                its(i) = r.Iterations(end);
            else
                its(i) = nan;
            end
        end
    end
    
end

function timing = plotLinearTimingsForExample(names, timing, its)
    bh = bar(timing, 'stacked');
%     bar(timing);
    % legend('Pre', 'Solve', 'Post', 'Total');
    legend('Solve', 'Pre', 'Post');
    set(gca, 'XTickLabel', names, 'TickLabelInterpreter', 'none');
    set(gca, 'XTickLabelRotation', -20)
    tot = sum(timing, 2);
    for i = 1:numel(tot)
        it = its(i);
        if it > 0
            s = num2str(it);
        elseif isnan(it)
            s = 'Fail';
        else
            s = '';
        end
        th = text(i, tot(i) + 0.05*max(tot), s);
        set(th, 'HorizontalAlignment', 'center', 'FontSize', 14)
    end
    ylim([0, 1.1*max(tot)]);
end

function d = getData(r)
    getter = @(x, fld) cellfun(@(x) x.(fld), x);
    prep = getter(r, 'PreparationTime');
    for i = 1:numel(r)
        x = r{i};
        if isfield(x, 'BlockAssembly')
            % Subtract block assembly, since this would normally be
            % measured as a part of the equations assembly (all other
            % solvers get a pre-assembled linear system in this benchmark)
            prep(i) = prep(i) - x.BlockAssembly;
        end
    end
    solve = getter(r, 'LinearSolutionTime');
    post = getter(r, 'PostProcessTime');
    tot = prep + solve + post;
    d = [solve; prep; post]';
end

%% Copyright Notice
%
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
