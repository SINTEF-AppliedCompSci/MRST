%% Large example demonstrating different types of fast assembly in MRST
% MRST has several different backends for automatic differentiation. In
% this example, we will show some of the different options for a large
% model with more than a million cells.
mrstModule add ad-core ad-blackoil ad-props mrst-gui coarsegrid linearsolvers
mrstModule add libgeometry deckformat

%% Build simulation model
% The simulation model is built manually. Here, we only give the necessary
% statements and refer to the original example in the MRST book for
% discussion and plots.
%
% By default, the model will only run a single step. To run the full
% schedule, set "singleStep = false" before running the script. Note that
% this will run three long schedules of 132 time-steps each, with a model
% with over a million cells!
if ~exist('singleStep', 'var')
    singleStep = true;
end
rng(0);
layers = [30 60 30];
dims = [100 100 sum(layers)];
[xmax,ymax, n]  = deal(1000*meter, 1000*meter, 30);
[x, y]  = meshgrid(linspace(0,xmax,n+1), linspace(0,ymax,n+1));
[x, y]  = deal(x',y');
dome    = 1-exp(sqrt((x - xmax/2).^2 + (y - ymax/2).^2)*1e-3);
[xn,yn] = deal(pi*x/xmax,pi*y/ymax);
perturb = sin(5*xn) + .5*sin(4*xn+6*yn) + cos(.25*xn*yn./pi^2) + cos(3*yn);
perturb = perturb/3.5;
[h, hr] = deal(8,1);
zt      = 50 + h*perturb + rand(size(x))*hr - 20*dome;
zb      = zt + 30;
zmb     = min(zb + 4 + 0.01*x - 0.020*y + hr*rand(size(x)), zb);
zmt     = max(zb -15 + 0.01*x - 0.025*y + hr*rand(size(x)), zt);

horizons = {struct('x', x, 'y', y, 'z', zt), ...
            struct('x', x, 'y', y, 'z', zmt), ...
            struct('x', x, 'y', y, 'z', zmb), ...
            struct('x', x, 'y', y, 'z', zb)};
grdecl   = convertHorizonsToGrid(horizons, 'dims', dims, 'layers', layers);

G = mprocessGRDECL(grdecl);
G = mcomputeGeometry(G);

% Petrophysics
rng(357371);
[K,L] = logNormLayers(G.cartDims, [100 400 10 50]*milli*darcy);
K     = K(G.cells.indexMap);
perm  = [K, K, 0.1*K];
rock  = makeRock(G, perm, 0.3);
%% Define wells and three-phase fluid model
% Define wells
simTime = 10*year;
pv      = poreVolume(G, rock);
injRate = 0.25*sum(pv)/simTime;
offset  = 10;

K = ceil(G.cartDims(3)*2/3);
nz = G.cartDims(3);
W = [];
W = verticalWell(W, G, rock, offset, offset, 1,...
                'Name', 'P1', 'comp_i', [1 0 0], 'sign', -1,...
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50);
W = verticalWell(W, G, rock,  offset, floor(G.cartDims(1)/2)+3, 1,...
                'Name', 'P2', 'comp_i', [1 0 0], 'sign', -1,...
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50);
W = verticalWell(W, G, rock, offset, G.cartDims(2) - offset/2, K, ...
                'Name', 'P3', 'comp_i', [1 0 0], 'sign', -1, ...
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50);
W = verticalWell(W, G, rock, G.cartDims(1)-5, offset, K,...
                'Name', 'I1', 'comp_i', [1 0 0], ...
                'Val', injRate, 'Type', 'rate', 'refDepth', 50);

% Three-phase template model with constant oil compressibility
fluid = initSimpleADIFluid('mu',    [1, 5, 0.1]*centi*poise, ...
                           'rho',   [1000, 700, 100]*kilogram/meter^3, ...
                           'c',     [1e-8, 1e-5, 1e-3]/barsa, ...
                           'n',     [2, 2, 1]);
p_ref    = 300*barsa;
% Construct reservoir model and initial state
gravity reset on
model = GenericBlackOilModel(G, rock, fluid);

region = getInitializationRegionsBlackOil(model, [90, 60]*meter, ...
            'datum_depth', 10*meter, 'datum_pressure', p_ref);
state0 = initStateBlackOilAD(model, region);

% Define simulation schedule and set solver parameters
timesteps = rampupTimesteps(simTime, 30*day, 10);
if singleStep
    timesteps = timesteps(1);
    casename = 'BlockAssemblyExampleSingleStep';
else
    casename = 'BlockAssemblyExample';
end
schedule = simpleSchedule(timesteps, 'W', W);

%% Plot initial setup
figure;
plotToolbar(G, state0);
plotWell(G, W);
%% Create packed problem handle
% The model is the only variation between the different solvers.
% Parametrize by the model.
packer = @(model, name, varargin) packSimulationProblem(state0, model, schedule, casename, 'name', name, varargin{:});
%% Sparse backend - default MRST setup
[model_sparse, nls] = setup_solver(model, 'sparse');
sparse_problem = packer(model_sparse, 'sparse-backend', 'NonLinearSolver', nls, 'description', 'Sparse');
%% Diagonal backend with MEX acceleration
% Faster version of AD that uses an intermediate diagonal representation
% and MEX acceleration to improve speed. Can generally replace the sparse
% version.
[model_diag, nls] = setup_solver(model, 'diagonal');
diag_problem = packer(model_diag, 'diagonal-backend', 'NonLinearSolver', nls, 'description', 'Diagonal');
%% Block-diagonal backend with MEX acceleration
% The fastest option, that assembles directly into AMGCL block matrices.
% Needs special linear solvers, and is limited in the complexity of
% non-reservoir equations. At the moment, the treatment in these
% solvers are only reliable when wells have a single cell each, or when the
% flow is driven by boundary conditions and sources.
[model_blockdiag, nls] = setup_solver(model, 'diagonal-block');
block_problem = packer(model_blockdiag, 'diagonal-backend-block', 'NonLinearSolver', nls, 'description', 'DiagonalBlock');
%% Simulate the problems
problems = {sparse_problem, diag_problem, block_problem};
simulatePackedProblem(problems);
%% Get output
[ws, states, reports, names] =...
    getMultiplePackedSimulatorOutputs(problems, 'readStatesFromDisk', false);
nd = numel(reports);
descr = cellfun(@(x) x.Description, problems, 'UniformOutput', false);
%% Plot timing per iteration/assembly
timings = cellfun(@(x) getReportTimings(x, 'total', true), reports);
d = arrayfun(@(x) [x.Assembly./x.NumberOfAssemblies, [x.LinearSolve, x.LinearSolvePrep]./x.Iterations],...
                    timings, 'UniformOutput', false);
figure;
bar(vertcat(d{:}))
set(gca, 'XTickLabel', names)
legend('Time per assembly', 'Time per linear solve', 'Prep per linear solve')
ylabel('Time [s]')
%% Plot as a function of different timesteps
t_per = cellfun(@(x) getReportTimings(x), reports, 'unif', false);
tot = cellfun(@(x) vertcat(x.Total), t_per, 'UniformOutput', false);
asm = cellfun(@(x) vertcat(x.Assembly) + vertcat(x.LinearSolvePrep), t_per, 'UniformOutput', false);
if numel(tot{1}) > 1
    figure('Position', [288, 695, 952, 283])
    hold on
    colors = lines(nd);
    tscale = 60;
    l = {};
    lw = 1.5;
    for i = 1:nd
        plot(cumsum(tot{i})/tscale, 'color', colors(i, :), 'linewidth', lw);
        l{end+1} = sprintf('%s: Total', names{i});
        plot(cumsum(asm{i})/tscale, '--', 'color', colors(i, :), 'linewidth', lw);
        l{end+1} = sprintf('%s: Assembly+Prep', names{i});
    end
    legend(l, 'location', 'northwest')
    axis tight
    xlabel('Timestep number')
    ylabel('Time [s]')
end
%% Setup routine used in script
function [model, nls] = setup_solver(model, type)
    nls = getNonLinearSolver(model);
    % We already have a short initial time-step
    nls.timeStepSelector.firstRampupStepRelative = 1;
    arg = {'tolerance', 1e-3, 'relaxation', 'ilu0', 'solver', 'bicgstab'};
    linsolve = AMGCL_CPRSolverAD(arg{:});
    linsolve.doApplyScalingCPR = false;
    linsolve.strategy = 'amgcl';
    nls.LinearSolver = linsolve;
    switch type
        case 'legacy'
            model = ThreePhaseBlackOilModel(model.G, model.rock, model.fluid, 'disgas', false, 'vapoil', false);
        case 'sparse'
            model.AutoDiffBackend = AutoDiffBackend();
        case 'diagonal'
            model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, 'rowMajor', true);
        case 'diagonal-block'
            model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, 'rowMajor', true, 'deferredAssembly', true);
            lsolve_block = AMGCL_CPRSolverBlockAD(arg{:});
            nls.LinearSolver = lsolve_block;
    end
    model = model.validateModel();
    ctm = model.FlowPropertyFunctions.ComponentTotalMass;
    ctm = ctm.setMinimumDerivatives([1e-10, 1e-6, 1e-6]);
    model.FlowPropertyFunctions.ComponentTotalMass = ctm;
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
