%% Large example demonstrating different types of fast assembly in MRST
% MRST has several different backends for automatic differentiation. In
% this example, we will show some of the different options for a large
% model with more than a million cells.
mrstModule add ad-core ad-blackoil ad-props mrst-gui coarsegrid linearsolvers
mrstModule add libgeometry deckformat
mrstModule add matlab_bgl

gravity off
%% Build simulation model
if ~exist('singleStep', 'var')
    singleStep = true;
end
if ~exist('dims', 'var')
    dims = [100, 100, 120];
end
if ~exist('mode', 'var')
    mode = 'lazy';
end
G = cartGrid(dims);
G = computeGeometry(G);
K     = 0.1*darcy;
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
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50, 'refDepth', -1);
W = verticalWell(W, G, rock,  offset, floor(G.cartDims(1)/2)+3, 1,...
                'Name', 'P2', 'comp_i', [1 0 0], 'sign', -1,...
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50, 'refDepth', -1);
W = verticalWell(W, G, rock, offset, G.cartDims(2) - offset/2, K, ...
                'Name', 'P3', 'comp_i', [1 0 0], 'sign', -1, ...
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50, 'refDepth', -1);
W = verticalWell(W, G, rock, G.cartDims(1)-5, offset, K,...
                'Name', 'I1', 'comp_i', [1 0 0], ...
                'Val', injRate, 'Type', 'rate', 'refDepth', 50, 'refDepth', -1);

% Three-phase template model with constant oil compressibility
fluid = initSimpleADIFluid('mu',    [1, 5, 0.1]*centi*poise, ...
                           'rho',   [1000, 700, 100]*kilogram/meter^3, ...
                           'c',     [1e-8, 1e-5, 1e-3]/barsa, ...
                           'n',     [2, 2, 1]);
p_ref    = 300*barsa;
% Construct reservoir model and initial state
gravity reset on
model = GenericBlackOilModel(G, rock, fluid);

state0 = initResSol(G, p_ref, [0, 0.5, 0.5]);

% Define simulation schedule and set solver parameters
timesteps = rampupTimesteps(simTime, 30*day, 10);
if singleStep
    timesteps = timesteps(1);
    caseName = 'BlockAssemblySimpleTestSingleStep';
else
    caseName = 'BlockAssemblySimpleTest';
end
schedule = simpleSchedule(timesteps, 'W', W);
%% Create packed problem handle
% The model is the only variation between the different solvers.
% Parametrize by the model.
packer = @(model, name, varargin) packSimulationProblem(state0, model, schedule, caseName, 'name', name, varargin{:});
%% Sparse backend - default MRST setup
[model_sparse, nls] = setupBlockTestSolver(model, 'sparse', mode);
sparse_problem = packer(model_sparse, 'sparse-backend', 'NonLinearSolver', nls, 'description', 'Sparse');
simulatePackedProblem(sparse_problem, 'restartStep', 1);

%% Diagonal backend with MEX acceleration
% Faster version of AD that uses an intermediate diagonal representation
% and MEX acceleration to improve speed. Can generally replace the sparse
% version.
[model_diag, nls] = setupBlockTestSolver(model, 'diagonal', mode);
diag_problem = packer(model_diag, 'diag-backend', 'NonLinearSolver', nls, 'description', 'Diagonal');
simulatePackedProblem(diag_problem, 'restartStep', 1);
%% Block-diagonal backend with MEX acceleration
% The fastest option, that assembles directly into AMGCL block matrices.
% Needs special linear solvers, and is limited in the complexity of
% non-reservoir equations. At the moment, the treatment in these
% solvers are only reliable when wells have a single cell each, or when the
% flow is driven by boundary conditions and sources.
[model_blockdiag, nls] = setupBlockTestSolver(model, 'diagonal-block', mode);
block_problem = packer(model_blockdiag, 'block-backend', 'NonLinearSolver', nls, 'description', 'DiagonalBlock');
simulatePackedProblem(block_problem, 'restartStep', 1);

%% Simulate the problems
problems = {sparse_problem, diag_problem, block_problem};
%% Get output
[ws, states, reports, names] =...
    getMultiplePackedSimulatorOutputs(problems, 'readStatesFromDisk', false);
nd = numel(reports);
descr = cellfun(@(x) x.Description, problems, 'UniformOutput', false);
%% Plot timing per iteration/assembly
timings = cellfun(@(x) getReportTimings(x, 'total', true), reports);
d = arrayfun(@(x) [x.Assembly./x.NumberOfAssemblies, [x.LinearSolve, x.LinearSolvePrep]./x.Iterations],...
                    timings, 'UniformOutput', false);
if ~mrstPlatform('gui')
    return
end
%% Plot time spent
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
