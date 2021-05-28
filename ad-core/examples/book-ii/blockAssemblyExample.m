%% Example demonstrating deferred assembly (WIP)
mrstModule add ad-core ad-blackoil
mrstModule add linearsolvers

gravity reset off
mrstModule add ad-core ad-props ad-blackoil
%%
if ~exist('forces', 'var')
    forces = 'pvwells';
end
if ~exist('dims', 'var')
    dims = [100, 100, 100];
end
if ~exist('singleStep', 'var')
    singleStep = true;
end

% Scaling for simple problem
muscale = centi*poise;
pscale = 100*barsa;
pdims = [1000, 1000, 100];
pvscale = 100000;
tscale = 60*day;
permscale = darcy;

G = cartGrid(dims, pdims);
G = computeGeometry(G);

rock = makeRock(G, permscale, 0.5);

fluid = initSimpleADIFluid('mu', [1, 1, 1]*muscale, 'rho', [100, 100, 100], 'cR', 1e-10/barsa, 'c', [1e-5, 1e-5, 1e-4]/barsa);
model = GenericBlackOilModel(G, rock, fluid, 'disgas', false, 'vapoil', false);
ncomp = 3;


isOil = G.cells.centroids(:, 3) > pdims(3)/3;

[so, sg] = deal(zeros(G.cells.num, 1));
p0 = pscale;
p = repmat(p0, G.cells.num, 1);

dp = 0.5*pscale;
dt = rampupTimesteps(364*5*tscale/60, tscale);
bc = [];
W = [];

e = 1e-3;
state0 = [];
schedule = [];
flux = sum(model.operators.pv)/sum(dt);
switch forces
    case 'bc'
        bc = fluxside(bc, G, 'XMax', flux, 'sat', [1-2*e, e, e]);
        bc = pside(bc, G, 'XMin', mean(p) - dp, 'sat', [1-2*e, e, e]);

    case 'wells'
        % ival = p0 + dp;
        % itype = 'bhp';
        ival = flux;
        itype = 'rate';
        W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', [1-2*e, e, e], 'type', itype, 'val', ival);
        W = verticalWell(W, G, rock, dims(1), dims(2), [], 'comp_i', [1-2*e, e, e], 'type', 'bhp', 'val', p0 - dp);
    case 'pvwells'
        [ii, jj, kk] = gridLogicalIndices(G);
        inj = ii == 1 & jj == 1;
        prod = ii == dims(1) & jj == dims(2);

        p(inj) = p(inj) + dp;
        p(prod) = p(prod) - dp;

        so(isOil) = 1 - e;
        sg(~isOil) = 1 - e;

        so(inj) = e;
        sg(inj) = e;
        model.operators.pv(inj | prod) = pvscale*model.operators.pv(inj | prod);
    case 'spe1'
        [G, rock, fluid, deck, state0] = setupSPE1();
        model = GenericBlackOilModel(G, rock, fluid, 'inputdata', deck);
        gravity reset on
        schedule = convertDeckScheduleToMRST(model, deck);
    otherwise
        error('Forces %s not supported', forces);
end

if isempty(state0)
    so = max(so, 1e-8);
    sg = max(sg, 1e-8);
    s0 = [1-so-sg, so, sg];
    state0 = initResSol(G, p, s0);
end
if singleStep
    dt = dt(1);
end
if isempty(schedule)
    schedule = simpleSchedule(dt, 'bc', bc, 'W', W);
end
% These are here so we can avoid running some of the cells below and still
% get plots.
[reportsparse, reportdiagcol, reportdiagrow, reportblock] = deal([]); %#ok

bz = ncomp;
use_cpr = true;
lsolve = pickLinearSolver(model, bz, use_cpr);
%% Sparse
modelsparse = model;
modelsparse.AutoDiffBackend = SparseAutoDiffBackend();
modelsparse = modelsparse.validateModel();

[~, statessparse, reportsparse] = simulateScheduleAD(state0, modelsparse, schedule, 'LinearSolver', lsolve);

%% Col major diag model
modeldiagcol = model;
modeldiagcol.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, 'rowMajor', false);
modeldiagcol = modeldiagcol.validateModel();

[~, statesdiagcol, reportdiagcol] = simulateScheduleAD(state0, modeldiagcol, schedule, 'LinearSolver', lsolve);
%% Row major diag model
modeldiagrow = model;
modeldiagrow.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, 'rowMajor', true);
modeldiagrow = modeldiagrow.validateModel();

[wsdiagrow, statesdiagrow, reportdiagrow] = simulateScheduleAD(state0, modeldiagrow, schedule, 'LinearSolver', lsolve);

%% Block model
modelblock = model;
modelblock.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, 'rowMajor', true);
modelblock.AutoDiffBackend.deferredAssembly = true;
modelblock = modelblock.validateModel();
modelblock.FacilityModel.toleranceWellRate = 1e-2;

lsolvenew = pickLinearSolver(modelblock, bz, use_cpr, true);
% lsolvenew.reductionStrategy = 'pre'
[wsblock, statesblock, reportblock] = simulateScheduleAD(state0, modelblock, schedule, 'LinearSolver', lsolvenew);
%% Plot results
reports = {reportsparse, reportdiagcol, reportdiagrow, reportblock};
names = {'Sparse (ADI)', 'DiagCol', 'DiagonalRow', 'DiagonalBlockRow'};

act = cellfun(@(x) ~isempty(x), reports);
reports = reports(act);
names = names(act);

results = cellfun(@(x) getReportTimings(x,  'skipConverged', true), reports, 'UniformOutput', false);
fn = @(x) [sum([x.Assembly]), sum([x.LinearSolvePrep]), sum([x.LinearSolve]), sum([x.Total])]./sum([x.Iterations]);
v = cellfun(fn, results, 'UniformOutput', false);
d = vertcat(v{:});
l = {'Equation assembly', 'Linear system assembly', 'Linear solve', 'Total time'};
figure(1); clf
% bar(d, 'stacked');
bar(d);

set(gca, 'XTickLabels', names, 'XTickLabelRotation', -15);
grid on
legend(l);
ylabel('Seconds per iteration')
t = sprintf('Grid %d x %d x %d (%d dof)', G.cartDims, G.cells.num*bz);
title(t)
%% Plot speedup
speedup = bsxfun(@rdivide, d(1, :), d);
clf;
bar(speedup(2:end, :))
legend(l);
ylabel(sprintf('Speedup relative to %s', names{1}))
set(gca, 'XTickLabels', names(2:end));

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

%%
function lsolve = pickLinearSolver(model, bz, use_cpr, new_solver)
    if nargin < 4
        new_solver = false;
    end
    b = model.AutoDiffBackend;
    if new_solver
        if use_cpr
            lsolve = AMGCL_CPRSolverBlockAD('strategy', 'amgcl');
        else
            lsolve = AMGCLSolverBlockAD();
        end
    else
        ord = getCellMajorReordering(model.G.cells.num, bz);
        if use_cpr
            lsolve = AMGCL_CPRSolverAD('block_size', bz, 'strategy', 'amgcl');
            lsolve.doApplyScalingCPR = false;
        else
            lsolve = AMGCLSolverAD('block_size', bz);
        end
        lsolve.equationOrdering = ord;
        lsolve.variableOrdering = ord;
    end
    
    if use_cpr
        lsolve.setSRelaxation('ilu0');
    else
        lsolve.setPreconditioner('relaxation');
        lsolve.setRelaxation('ilu0');
    end
    lsolve.tolerance = 1e-3;
%     lsolve.verbose = true;
%     lsolve.amgcl_setup.verbose = 100;
end
