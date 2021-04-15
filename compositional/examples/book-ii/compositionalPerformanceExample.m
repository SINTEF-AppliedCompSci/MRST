%% Conceptual example on how to run larger compositional problems
% The default settings of the compositional module are suitable for small
% problems and should work on most configurations with a recent MATLAB
% version. However, with an external linear solver and a change of backend,
% larger problems can be solved.
%
% The example is discussed in Section 8.5.5 of the second MRST book:
% Advanced Modelling with the MATLAB Reservoir Simulation Toolbox (MRST),
% Cambridge University Press, 2021.
mrstModule add compositional ad-core linearsolvers ad-props ad-blackoil

%% Set up the model problem used as benchmark
% The default benchmark is specified on a 50x50x50 Cartesian grid with a
% single producer-injector pair and a six-component fluid model taken from
% the 5th SPE Comparative Solution Project, which gives a total of 750 000
% reservoir degrees of freedom. This system will be used to compare the
% computational efficiency of various backends. If you want another test
% case, it is straightforward to change the fluid mixture or modify the
% grid and then perform the test on your own computer.

% Grid and petrophysics
gravity reset on; rng(0)
[nx,ny,nz] = deal(50);
dims  = [nx, ny, nz];
pdims = [1000, 1000, 100];
G     = computeGeometry(cartGrid(dims, pdims));
K     = logNormLayers(G.cartDims);
rock  = makeRock(G, K*milli*darcy, 0.2);
pv    = poreVolume(G, rock);

% Fluid model
[fluid, info] = getBenchmarkMixture('spe5');
flowfluid = initSimpleADIFluid('rho', [1000, 500, 500], ...
                       'mu', [1, 1, 1]*centi*poise, ...
                       'n', [2, 2, 2], ...
                       'c', [1e-5, 0, 0]/barsa);

% Drive mechanisms
minP    = 0.5*info.pressure;
totTime = 10*year;
W = verticalWell([], G, rock, 1, 1, 1, 'comp_i', [1, 0], ...
    'Type', 'rate', 'name', 'Inj', 'Val', sum(pv)/totTime, 'sign', 1);
W = verticalWell(W, G, rock, nx, ny, nz, ...
    'comp_i', [0.5, 0.5], 'Name', 'Prod', 'Val', minP, 'sign', -1);
for i = 1:numel(W)
    W(i).components = info.injection;
end
bc = [];    % No boundary conditions

%% Initialize models
% We initialize four different models: The first uses the standard
% constructor with the Sparse backend for AD. The second uses the Diagonal
% AD backend instead, which is faster for problems with a moderate to many
% degrees of freedom and several components. The third extends the diagonal
% backend with C++-acceleration through a mex interface. A compiler must be
% available (see mex -setup) for this to work. Note that the first
% simulation with diagonal+mex will occasionally stop to compile
% subroutines, so the user is encouraged to repeat the benchmark if this is
% the case. The fourth stores the nonzero diagonals of the Jacobi matrices
% using row-major rather than MATLAB's default column-major format. It also
% uses deferred assembly. (See Chapter 6 in the second MRST book for more
% details.)
%
% The performance of the backends varies with your computational
% configuration, so we generally encourage you to test different versions
% to determine which option gives the best computational performance for
% your model and computer.
arg = {G, rock, flowfluid, fluid, 'water', false};
diagonal_backend = DiagonalAutoDiffBackend();
mex_backend      = DiagonalAutoDiffBackend('useMex', true);
row_backend      = DiagonalAutoDiffBackend('useMex', true, 'rowMajor', true, 'deferredAssembly', true);
sparse_backend   = SparseAutoDiffBackend();

constructor = @GenericOverallCompositionModel;

modelSparseAD      = constructor(arg{:}, 'AutoDiffBackend', sparse_backend);
modelDiagonalAD    = constructor(arg{:}, 'AutoDiffBackend', diagonal_backend);
modelMexDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', mex_backend);
modelRowDiagonalAD = constructor(arg{:}, 'AutoDiffBackend', row_backend);

%% Set up driving forces and initial state
ncomp  = fluid.getNumberOfComponents();
s0     = [1, 0];
state0 = initCompositionalState(G, info.pressure, info.temp, s0, ...
                                info.initial, modelSparseAD.EOSModel);

%% Pick linear solver
% The AMGCL library is one possible solver option for MRST. It is fairly
% easy to write interfaces to other solvers using MEX files and/or the
% LinearSolverAD class. We call the routine for automatically selecting a
% reasonable linear solver for the specific model.
linsolve = selectLinearSolverAD(modelDiagonalAD);
disp(linsolve)

%% Solve a single, short step to benchmark assembly and solve time
shortSchedule  = simpleSchedule(0.1*day, 'bc', bc, 'W', W);
nls            = NonLinearSolver('LinearSolver', linsolve);
linsolve_block = AMGCL_CPRSolverBlockAD('tolerance', 1e-4);

[~, statesSparse, reportSparse] = ...
    simulateScheduleAD(state0, modelSparseAD, shortSchedule, 'nonlinearsolver', nls);
[~, statesDiagonal, reportDiagonal] = ...
    simulateScheduleAD(state0, modelDiagonalAD, shortSchedule, 'nonlinearsolver', nls);
[~, statesMexDiagonal, reportMexDiagonal] = ...
    simulateScheduleAD(state0, modelMexDiagonalAD, shortSchedule, 'nonlinearsolver', nls);
nls.LinearSolver = linsolve_block;
[~, ~, reportRowDiagonal] = ...
    simulateScheduleAD(state0, modelRowDiagonalAD, shortSchedule, 'nonlinearsolver', nls);

%% Plot the CPU-time statistics 
figure(1); clf
time_sparse   = getReportTimings(reportSparse);
time_diagonal = getReportTimings(reportDiagonal);
time_row      = getReportTimings(reportRowDiagonal);
time_mex      = getReportTimings(reportMexDiagonal);
get = @(x) [x.Assembly./x.NumberOfAssemblies, ...
            (x.LinearSolve + x.LinearSolvePrep)./x.Iterations, ...
            x.Total./x.Iterations];
time = [get(time_sparse); get(time_diagonal); get(time_mex); get(time_row)];
bar(time)
set(gca, 'XTickLabel', {'Sparse', 'Diagonal', 'Diagonal (MEX)', ...
    'Diagonal (MEX+RowMajor)'},'XTickLabelRotation',10);
legend('Assembly / each', 'Linear solve / it', 'Total / it', 'Location', 'NorthEast')
ylabel('Time [s]')

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
