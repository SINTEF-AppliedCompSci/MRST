%% Example demonstrating accelerated assembly for faster simulation
% Running larger cases in MRST is possible, but the default parameters of
% MRST's solvers are optimized for ease-of-use out of the box with just a
% standard Matlab installation. We do not wish to require users to have C++
% compilers, external software or expensive add-on toolboxes to make use of
% the software.
%
% However, simulation can be achieved much faster if one has access to both
%  - A C++ compiler
%  - An external linear solver with a Matlab interface.
%
% This examples demonstrates the use of MEX to seamlessly improve assembly
% speed for AD-solvers.
mrstModule add ad-core ad-blackoil deckformat ad-props linearsolvers
%% Define a 70x70x30 grid with uniform properties
% Obviously, you can adjust this to test your setup
dims = [70, 70, 30];
pdims = [1000, 1000, 100]*meter;
G = cartGrid(dims, pdims);
G = computeGeometry(G);
rock = makeRock(G, 0.1*darcy, 0.3);
%% Set up QFS water injection scenario
% As the point of this example is on how to set up the solvers, we just set
% up a simple injection scenario for a 3D quarter-five-spot problem.
gravity reset on
time = 10*year;
irate = 0.3*sum(poreVolume(G, rock)/time);

W = [];
W = verticalWell(W, G, rock, 1, 1, [], 'comp_i', [1, 0, 0], ...
                'type', 'rate', 'val', irate);
W = verticalWell(W, G, rock, dims(1), dims(2), 1, 'comp_i', [1, 0, 0], ...
                 'type', 'bhp', 'val', 100*barsa);
%% Three-phase fluid model with compressibility and nonlinear flux
fluid = initSimpleADIFluid('n', [2, 3, 1], ...
                           'rho', [1000, 700, 10], ...
                           'mu', [1, 5, 0.1]*centi*poise, ...
                           'c', [0, 1e-6, 1e-4]/barsa);

%% Do 5 time-steps of 30 days each to test assembly speed
dt = repmat(30*day, 5, 1);
% Alternatively: Run the whole schedule with 200 time-steps.
% dt = rampupTimesteps(time, time/200);
schedule = simpleSchedule(dt, 'W', W);

model = GenericBlackOilModel(G, rock, fluid, 'disgas', false, 'vapoil', false);
%% Set up mex-accelerated backend and reduced variable set for wells
% We use the Diagonal autodiff backend to calculate derivatives. By
% default, this uses only Matlab, so we also set the "useMex" flag to be
% true. It will then use C++ versions of most discrete operators during
% assembly.
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
model = model.validateModel();
model.FacilityModel.primaryVariableSet = 'bhp';
%% Run assembly tests (will compile backends if required)
% Note: Benchmarks are obviously not accurate if compilation is performed.
testMexDiagonalOperators(model, 'block_size', 3);
%% Set up a compiled linear solver
% We use a AMGCL-based CPR solver to solve the problem. It is sufficient to
% have a working C++ compiler, the AMGCL repository and BOOST available to
% compile it. See the documentation of amgcl_matlab for more details on how
% to set up these paths.

ncomp = model.getNumberOfComponents();
solver = AMGCL_CPRSolverAD('tolerance', 1e-3, 'block_size', ncomp, ...
                           'coarsening', 'aggregation', 'relaxation', 'ilu0');

nls = NonLinearSolver('LinearSolver', solver);
%% Set up initial state
% We define fluid contacts and datum pressure
region = getInitializationRegionsBlackOil(model, [70, 30], 'datum_pressure', 200*barsa);
state0 = initStateBlackOilAD(model, region);

figure(1); clf
plotCellData(G, state0.s(:, [2, 3, 1]));
view(30, 30);
plotWell(G, W);
title('Initial saturations');
%% Run a serial simulation
maxNumCompThreads(1)
solver.amgcl_setup.nthreads = 1; % Specify single-thread manually
[ws0, states0, rep0] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);
%% Run a parallel simulation with four threads
maxNumCompThreads(4)
solver.amgcl_setup.nthreads = 4; % Specify threads manually
[ws, states, rep] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);
% Reset max number of threads
N = maxNumCompThreads('automatic');
%% Plot the final saturations
figure(1); clf
plotCellData(G, states{end}.s(:, [2, 3, 1]));
view(30, 30);
plotWell(G, W);
title('Final saturations');
%% Plot timings
% We plot the total time, the assembly (linearization) time, the linear
% solver time and the linear solver preparation time. The assembly and
% linear solver time are in general parallel, while the preparation is
% currently not.
%
% We likely observe that the linear solver is the place where most of the
% time is spent: Manipulating sparse matrices to fit them to an external
% solver can take significant time for a moderate-size problem.
serial = getReportTimings(rep0);
par = getReportTimings(rep);

total0 = sum(vertcat(serial.Total));
assembly0 = sum(vertcat(serial.Assembly));
lsolve0 = sum(vertcat(serial.LinearSolve));
lsolve_prep0 = sum(vertcat(serial.LinearSolvePrep));
its0 = rep0.Iterations;

total = sum(vertcat(par.Total));
assembly = sum(vertcat(par.Assembly));
lsolve = sum(vertcat(par.LinearSolve));
lsolve_prep = sum(vertcat(par.LinearSolvePrep));
its = rep.Iterations;

a = [total0, assembly0, lsolve0, lsolve_prep0];
b = [total, assembly, lsolve, lsolve_prep];

a(end+1) = a(1) - sum(a(2:end));
b(end+1) = b(1) - sum(b(2:end));

figure;
bar([a; b])
legend('Total time', 'Assembly', 'Linear solve', 'Linear solve preparation', 'Other functions')
set(gca, 'XTickLabel', {'Serial', [num2str(N), ' threads']})
%% Other notes
% We note that the more degrees-of-freedom you have per cell in your model,
% the more efficient the MEX operators are compared to Matlab. See e.g. the
% compositional module's example "compositionalLargerProblemTutorial.m"
% for another example benchmark.

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
