%% Conceptual example on how to run larger compositional problems
% The default settings of the compositional module are suitable for small
% problems and should work on most configurations with a recent Matlab
% version. However, with an external linear solver and a change of backend,
% larger problems can be solved.
mrstModule add compositional ad-core linearsolvers ad-props
useBC = false; % Use BC instead of wells
includeWater = false; % Include aqueous phase
if ~exist('useNatural', 'var')
    useNatural = true; % Use natural variables formulation
end
% Define a problem
gravity reset on
nx = 30;
ny = 30;
nz = 30;

dims = [nx, ny, nz];
pdims = [1000, 1000, 100];
G = cartGrid(dims, pdims);
G = computeGeometry(G);
% Random permeability field
rng(0);
K = logNormLayers(G.cartDims);
K = K(G.cells.indexMap);

rock = makeRock(G, K*milli*darcy, 0.2);
pv = poreVolume(G, rock);
% Take the SPE5 fluid model (six components)
[fluid, info] = getBenchmarkMixture('spe5');

minP = 0.5*info.pressure;
totTime = 10*year;

% Set up driving forces
[bc, W] = deal([]);
if useBC
    bc = fluxside(bc, G, 'xmin', sum(pv)/totTime, 'sat', [0, 1, 0]);
    bc = pside(bc, G, 'xmax', minP, 'sat', [0, 1, 0]);

    bc.components = repmat(info.injection, numel(bc.face), 1);
else
    W = verticalWell(W, G, rock, 1, 1, 1, 'comp_i', [0, 1, 0], ...
        'Type', 'rate', 'name', 'Inj', 'Val', sum(pv)/totTime, 'sign', 1);

    W = verticalWell(W, G, rock, nx, ny, nz, ...
        'comp_i', [0, 0.5, 0.5], 'Name', 'Prod', 'Val', minP, 'sign', -1);


    for i = 1:numel(W)
        W(i).components = info.injection;
    end
end
flowfluid = initSimpleADIFluid('rho', [1000, 500, 500], ...
                       'mu', [1, 1, 1]*centi*poise, ...
                       'n', [2, 2, 2], ...
                       'c', [1e-5, 0, 0]/barsa);
%% Initialize models
% We initialize two models: The first uses the standard constructor, and
% uses the Sparse backend for AD. The second example uses the Diagonal AD
% backend instead, which is faster for problems with a moderate to many
% degrees of freedom and several components. The third option is to
% extended the diagonal backend with C++-acceleration through a mex
% interface. A compiler must be available (see mex -setup) for this to
% work. Note that the first simulation with diagonal+mex will occasionally
% stop to compile subroutines, so the user is encouraged to repeat the
% benchmark if this is the case.
%
% The performance of the backends varies from configuration to
% configuration, so the user is encouraged to test different versions until
% the best speed is found.
arg = {G, rock, flowfluid, fluid, 'water', includeWater};
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true, 'rowMajor', true);
sparse_backend = SparseAutoDiffBackend();

if useNatural
    constructor = @GenericNaturalVariablesModel;
else
    constructor = @GenericOverallCompositionModel;
end

modelSparseAD = constructor(arg{:}, 'AutoDiffBackend', sparse_backend);
modelDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', diagonal_backend);
modelMexDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', mex_backend);

%% Set up driving forces and initial state
ncomp = fluid.getNumberOfComponents();
if includeWater
    s0 = [0.2, 0.8, 0];
else
    s0 = [1, 0];
    for i = 1:numel(W)
        W(i).compi = W(i).compi(2:end);
    end
    if ~isempty(bc)
        bc.sat = bc.sat(:, 2:end);
    end
end
state0 = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, modelSparseAD.EOSModel);
%% Pick linear solver
% The AMGCL library is one possible solver option for MRST. It is fairly
% easy to write interfaces to other solvers using MEX files and/or the
% LinearSolverAD class. We call the routine for automatically selecting a
% reasonable linear solver for the specific model.
linsolve = selectLinearSolverAD(modelDiagonalAD);
disp(linsolve)
%% Solve a single short step to benchmark assembly and solve time
shortSchedule = simpleSchedule(0.1*day, 'bc', bc, 'W', W);

nls = NonLinearSolver('LinearSolver', linsolve);
[~, statesSparse, reportSparse] = simulateScheduleAD(state0, modelSparseAD, shortSchedule, 'nonlinearsolver', nls);
[~, statesDiagonal, reportDiagonal] = simulateScheduleAD(state0, modelDiagonalAD, shortSchedule, 'nonlinearsolver', nls);
try
    [~, statesMexDiagonal, reportMexDiagonal] = simulateScheduleAD(state0, modelMexDiagonalAD, shortSchedule, 'nonlinearsolver', nls);
catch
    [statesMexDiagonal, reportMexDiagonal] = deal([]);
end

%% Solve using direct solver
% This system is too large for the standard direct solver in Matlab. This
% may take some time!
[~, ~, reportDirect] = simulateScheduleAD(state0, modelSparseAD, shortSchedule);

%% Plot the time taken to solve a single step
figure(1 + useNatural); clf

getTime = @(report) [sum(cellfun(@(x) x.AssemblyTime, report.ControlstepReports{1}.StepReports{1}.NonlinearReport)), ... % Assembly
                     sum(cellfun(@(x) x.LinearSolver.SolverTime, report.ControlstepReports{1}.StepReports{1}.NonlinearReport(1:end-1))),... % Linear solver
                     sum(report.SimulationTime), ...
                     ];
time_sparse = getTime(reportSparse);
time_diagonal = getTime(reportDiagonal);
time_direct = getTime(reportDirect);
if isempty(reportMexDiagonal)
    time = [time_sparse; time_diagonal; time_direct];
    bar(time)
    set(gca, 'XTickLabel', {'Sparse', 'Diagonal', 'Sparse+Direct solver'});
else
    time_mex = getTime(reportMexDiagonal);
    time = [time_sparse; time_diagonal; time_mex; time_direct];
    bar(time)
    set(gca, 'XTickLabel', {'Sparse', 'Diagonal', 'Diagonal with Mex', 'Sparse+Direct solver'});
end
legend('Equation assembly', 'Linear solver', 'Total time', 'Location', 'NorthWest')
%% Zoom in on assembly time
assembly_time = time(:, 1);
ylim([0, 1.2*max(assembly_time)]);

%%
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
