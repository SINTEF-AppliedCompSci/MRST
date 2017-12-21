%% Conceptual example on how to run larger compositional problems
% The default settings of the compositional module are suitable for small
% problems and should work on most configurations with a recent Matlab
% version. However, with an external linear solver and a change of backend,
% larger problems can be solved.
mrstModule add compositional ad-core linearsolvers ad-props
useBC = false; % Use BC instead of wells
includeWater = false; % Include aqueous phase

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
[fluid, info] = getCompositionalFluidCase('spe5');

minP = 0.5*info.pressure;
totTime = 10*year;

% Set up driving forces
[bc, W] = deal([]);
if useBC
    bc = fluxside(bc, G, 'xmin', sum(pv)/totTime, 'sat', [0, 1, 0]);
    bc = pside(bc, G, 'xmax', minP, 'sat', [0, 1, 0]);

    bc.components = repmat(info.injection, numel(bc.face), 1);
else
    W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', [0, 1, 0], ...
        'Type', 'rate', 'name', 'Inj', 'Val', sum(pv)/totTime, 'sign', 1);

    W = verticalWell(W, G, rock, nx, ny, [], ...
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
% degrees of freedom and several components.
%
% The performance of the backends varies from configuration to
% configuration, so the user is encouraged to test different versions until
% the best speed is found.
arg = {G, rock, flowfluid, fluid, 'water', includeWater};
modelSparseAD = NaturalVariablesCompositionalModel(arg{:});
modelDiagonalAD = NaturalVariablesCompositionalModel(arg{:}, 'AutoDiffBackend', DiagonalAutoDiffBackend('modifyOperators', true));

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
% LinearSolverAD class.
useAMGCL = true;
if useAMGCL
    mrstModule add linearsolvers
    linsolve = AMGCLSolverAD('preconditioner', 'relaxation', ...
                             'relaxation', 'ilu0', ...
                             'solver', 'bicgstab', ...
                             'maxIterations', 100, ...
                             'tolerance', 1e-3);
else
    linsolve = GMRES_ILUSolverAD('tolerance', 1e-3, 'maxIterations', 100);
end
% The natural variables solver reduces linear systems internally, so we
% cannot use the reduceToCell option.
linsolve.reduceToCell = false;
% Instead, we only keep the cell variables using the low-level interface.
linsolve.keepNumber = (ncomp+modelDiagonalAD.water)*G.cells.num;
% Right diagonal scaling helps with the condition number of the system due
% to mixed derivatives
linsolve.applyRightDiagonalScaling = true;

%% Solve a single short step to benchmark assembly and solve time
shortSchedule = simpleSchedule(0.1*day, 'bc', bc, 'W', W);

nls = NonLinearSolver('LinearSolver', linsolve);
[~, statesSparse, reportSparse] = simulateScheduleAD(state0, modelSparseAD, shortSchedule, 'nonlinearsolver', nls);
[~, statesDiagonal, reportDiagonal] = simulateScheduleAD(state0, modelDiagonalAD, shortSchedule, 'nonlinearsolver', nls);

%% Solve using direct solver
% This system is too large for the standard direct solver in Matlab. This
% may take some time!
[~, ~, reportDirect] = simulateScheduleAD(state0, modelSparseAD, shortSchedule);

%% Plot the time taken to solve a single step
figure(1); clf

getTime = @(report) [sum(cellfun(@(x) x.AssemblyTime, report.ControlstepReports{1}.StepReports{1}.NonlinearReport)), ... % Assembly
                     sum(cellfun(@(x) x.LinearSolver.SolverTime, report.ControlstepReports{1}.StepReports{1}.NonlinearReport(1:end-1))),... % Linear solver
                     sum(report.SimulationTime), ...
                     ];
time_sparse = getTime(reportSparse);
time_diagonal = getTime(reportDiagonal);
time_direct = getTime(reportDirect);


bar([time_sparse; time_diagonal; time_direct])
set(gca, 'XTickLabel', {'Sparse', 'Diagonal', 'Sparse+Direct solver'});
legend('Equation assembly', 'Linear solver', 'Total time', 'Location', 'NorthWest')

