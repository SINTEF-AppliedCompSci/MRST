%% Example demonstrating AMGCL on a few test problems
mrstModule add msrsb incomp coarsegrid spe10 ...
   linearsolvers ad-core ad-blackoil ad-props example-suite

%% Setup problem
testcase = 'spe10';
switch lower(testcase)
    case 'simple'
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
    case 'spe10'
        mrstModule add spe10
        [state0, model, schedule]  = setupSPE10_AD('layers', 1);
        forces = schedule.control(1);
        dt = schedule.step.val(1);
        G = model.G;
    case 'spe9'
        [G, rock, fluid, deck, state0] = setupSPE9();
        model = selectModelFromDeck(G, rock, fluid, deck);
        schedule = convertDeckScheduleToMRST(model, deck);
        forces = schedule.control(1);
        dt = schedule.step.val(1);
end
% Set up model to get a linearized system
model = model.validateModel(forces);
state0 = model.validateState(state0);
% Get equations for time-step
state = state0;
problem = model.getEquations(state0, state, dt, forces, 'iteration', inf);
[A0, b0] = problem.getLinearSystem();
% Use linear solver class to perform schur complement and remove wells from
% system
ncomp = model.water + model.oil + model.gas;
lsolve = BackslashSolverAD();
lsolve.keepNumber = ncomp*G.cells.num;
[A0, b0, sys] = lsolve.reduceLinearSystem(A0, b0);

% Make a simple pressure equation via an unweighted sum of the equations.
% In general, this is not always accurate. The AMGCL_CPRSolverAD solver
% class has implemented more sophisticated variants (e.g. dynamic row-sum,
% true-impes, quasi-impes, etc) which are used later in this example.
pix = 1:G.cells.num;
for i = 2:ncomp
    ix = (1:G.cells.num) + (i-1)*G.cells.num;
    A0(pix, :) = A0(pix, :) + A0(ix, :);
    b0(pix) = b0(pix) + b0(ix);
end
% Reorder to cell-major format (required by solver). MRST uses
% equation-major ordering by default.
subs = getCellMajorReordering(G.cells.num, ncomp);

b0 = b0./norm(b0, inf);
A = A0(subs, subs);
b = b0(subs);
% b = b0;
its = 100;
% AMGCL/C++ uses cell-major ordering, we do this outside to get more
% reliable timing results.
At = A';

mrstVerbose on
%% Call solver
tic();
x1 = callAMGCL_cpr(At, b, ncomp, 'isTransposed', true, 'maxIterations', its, 'cellMajorOrder', true);
t_wrapper = toc();  
%% Setup SPE9 black-oil case
[G, rock, fluid, deck, state] = setupSPE9();
model = selectModelFromDeck(G, rock, fluid, deck);
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
schedule = convertDeckScheduleToMRST(model, deck);
%% Solve the schedule with AMGCL-CPR as the linear solver
lsolve = AMGCL_CPRSolverAD('maxIterations', 200, 'tolerance', 1e-3);
lsolve.setSRelaxation('ilu0');
% We can set coarsening and solver options as well
lsolve.setCoarsening('aggregation');
lsolve.setSolver('bicgstab');
[wellSols, states, report] = simulateScheduleAD(state, model, schedule, 'linearsolver', lsolve);
%% Plot breakdown of simulator time vs linear solver time
nstep = numel(report.ControlstepReports);
ls_time = zeros(nstep, 1);
for i = 1:nstep
    rr = report.ControlstepReports{i};
    for j = 1:numel(rr.StepReports{1}.NonlinearReport)
        rrr = rr.StepReports{1}.NonlinearReport{j};
        if rrr.Converged
            continue
        end
        if isfield(rrr.LinearSolver, 'SolverTime')
            ls_time(i) = ls_time(i) + rrr.LinearSolver.SolverTime;
        end
    end
end
% Plot per-step breakdown
figure;
bar([ls_time, report.SimulationTime - ls_time], 'stacked')
legend('Linear solver time', 'Assembly');
xlabel('Timestep number');
ylabel('Time [s]');

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
