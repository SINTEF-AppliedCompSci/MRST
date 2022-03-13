%% Introduction to the use of multiscale solvers as an iterative two-level method
% This example demonstrates the use of multiscale solvers as two-level
% multigrid-like methods. We demonstrate several different interfaces for
% iterative multiscale methods in MRST, depending on the intended usage. We
% recommend that you first go through the introMultiscale example before
% running this example, as we here do not detail the basic interfaces.
mrstModule add coarsegrid spe10 msrsb incomp ad-core

%% Define grid and rock
dims = [30, 110];
physDims = dims .*[20, 10]*ft;
% Create grid
G = cartGrid(dims, physDims);
G = computeGeometry(G);
n = G.cells.num;

% Take the rock structure from layer 15
rock = getSPE10rock(1:dims(1), 1:dims(2), 15);
rock.perm = rock.perm(:, 1);
rock.poro = max(rock.poro, 1e-3);
hT = computeTrans(G, rock);
% Coarsen into uniform grid
cdims = ceil(dims./[5, 10]);
p = partitionUI(G, cdims);
% Make coarse grid
CG = generateCoarseGrid(G, p);
% Add coarse geometry (centroids, normals, etc)
CG = coarsenGeometry(CG);
CG = storeInteractionRegion(CG);

%% Set up a problem with two wells
% We here set up a problem where wells are the driving force. Wells are
% the most common driving forces for field-scale simulation, and can be
% challenging to resolve as they represent a highly localized source terms.
% Here, the wells are simple point-wells in the first and last cell of the
% domain, corresponding to each corner of a structured grid.
W = [];
W = addWell(W, G, rock, 1, 'type', 'bhp', 'val', 100*barsa);
W = addWell(W, G, rock, n, 'type', 'bhp', 'val', 50*barsa);
fluid = initSingleFluid('rho', 1, 'mu', 1);
state0 = initResSol(G, 0);
% Call fine-scale solver, with matrix output in state enabled
state = incompTPFA(state0, G, hT, fluid, 'MatrixOutput', true, 'W', W);

%% Get the linear system from the report
% We perform a Schur complement to remove the well equations. We are doing
% this since we are going to use the low-level iterative solver interface,
% where it is expected that the linear system is n by n in size, where n is
% the number of cells in the domain.
A = state.A;
b = state.rhs;
[A, b, A_ww, A_wp, q_w] = eliminateWellEquations(A, b, n);
% To recover the extra equations, for some solution p, you can use
% p = recoverWellSolution(A_ww, A_wp, q_w, p)

%% Create multiscale basis
% Build MsRSB basis functions
basis = getMultiscaleBasis(CG, A);

%% Set solver parameters
% Solver handle for the coarse scale solver
coarse_solver = @(A_c, b_c) mldivide(A_c, b_c);
% Initial guess
p0 = [];
% Verbosity enabled
verb = true;
% Tolerance
tol = 1e-8;
% Set the maximum number of cycles
maxIterations = 100;

%% Get second-stage relaxation
% We select a single cycle of ILU(0) as our second stage smoother. This
% routine provides us a function interface @(A, b) which itself returns a
% function handle. This final function handle is then an approximate
% inverse of A. We use this interface since it allows for a setup phase,
% for example to perform a partial factorization of A. ILU(0) is one such
% partial factorization, where an incomplete lower-upper factorization is
% performed, allowing L and U to not have more non-zero elements than A.
fn = getSmootherFunction('type', 'ilu0', 'iterations', 1);
% Perform ILU0 factorization (using Matlab's builtin routines)
prec = fn(A, b);
p_approx = prec(b);

%% Plot the action of the preconditioner
% We plot the result of a single application of ILU(0). We observe that the
% preconditioner only modifies the solution very close to the wells. By
% inheriting the sparsity pattern of the system matrix, which is a local
% discretization where spatial neighbors depend only on each other, a
% single application of the preconditioner will not account for the global
% nature of the elliptic pressure equation. In the parlance of multigrid
% and multilevel methods, ILU(0) is a smoother, as it is quite efficient at
% removing local (or high-frequency) errors, but not quite accurate when it
% comes to global (or low-frequency) errors.
clf;
subplot(1, 2, 2)
plotCellData(G, state.pressure, 'EdgeColor', 'none');
axis equal tight, cax=caxis();
title('Solution');
subplot(1, 2, 1)
plotCellData(G, p_approx, 'EdgeColor', 'none');
axis equal tight, caxis(cax);
title('One pass of ILU(0)');

%% Solve and plot a regular multiscale solution
% We compute a simple non-iterative multiscale approximation of the
% pressure field. We plot the solution together with the reference.
state_ms = incompMultiscale(state, CG, hT, fluid, basis, 'W', W);
subplot(1, 2, 1), cla
plotCellData(G, state_ms.pressure, 'EdgeColor', 'none');
axis equal tight, caxis(cax);
title('Regular multiscale solution');

%% Plot residual error in the multiscale solution
% We plot the cell-wise absolute residual error,
% $$ |A \mathbf{p} - b| $$
% for the multiscale solver. We observe that the error is highly localized
% around the support of each basis function, forming what is essentially
% the contours of a dual grid. This residual error is highly localized and
% is the result of the limited support of the basis functions, which are
% imposed to ensure that the coarse-scale system is sparse in nature.
subplot(1, 2, 2); cla
plotCellData(G, abs(A*state_ms.pressure - b), 'EdgeColor', 'none');
plotGrid(CG, 'EdgeColor', 'w', 'linewidth', 1, 'FaceColor', 'none')
caxis([0, 1e-10])
axis equal tight
title('Residual error, MS');

%% Perform two solutions: One without GMRES and one with GMRES.
% The interplay between local smoothers and coarse-grid corrections is
% widely studied for multi-level solvers. The central idea is to use the
% inexpensive smoother to remove local errors, and the coarse scale solver
% to remove low-frequency errors. Together, the two solvers can fairly
% efficiently remove error modes from elliptic or near-elliptic problems.
%
% We use the low-level solveMultiscaleIteratively interface here.
% When useGMRES is false, a simple preconditioned Richardson iteration is
% performed. This problem is highly amenable to Krylov acceleration, since
% the error modes span multiple scales.
useGMRES = false;
[p_ms, report] = solveMultiscaleIteratively(A, b, p0, basis, fn, tol, ....
    maxIterations, coarse_solver, useGMRES, verb);
% Solve with GMRES
useGMRES = true;
[~, report_gmres] = solveMultiscaleIteratively(A, b, p0, basis, fn, tol, ...
    maxIterations, coarse_solver, useGMRES, verb);

%% Plot the residual evolution throughout the iterative solver
% We clearly observe that GMRES is much more efficient than the simple
% Richardson solver when measured in the number of iterations.
clf; hold on;
plot(report.resvec,'-or','MarkerFaceColor',[.8 .8 .8]);
plot(report_gmres.resvec,'-ob','MarkerFaceColor',[.8 .8 .8]);
set(gca, 'YScale', 'log')
legend('MS+ILU(0)', 'MS+ILU0, GMRES');
xlabel('Iterations');
ylabel('Residual')

%% We can also achieve the same result through the incompMultiscale interface
% This is the easiest way of improving a multiscale solution. Here, we
% solve only five iterations in order to inexpensively improve upon the
% initial approximation.

state_ms_it = incompMultiscale(state, CG, hT, fluid, basis, 'W', W, ...
    'getSmoother', fn, 'iterations', 5, 'useGMRES', true, 'tolerance', tol);

%% Plot discrete divergence of iterated solution
% The goal of iterative multiscale methods is not necessarily to solve the
% solution to a very strict tolerance, but rather to get an acceptable
% pressure solution together with a divergence-free velocity field. Here,
% we plot the divergence of both the fine-scale and the iterated multiscale
% solutions to see that both are non-zero only in the well cells, even when
% the solution is not fully converged.
op = setupOperatorsTPFA(G, rock);

Div = @(flux) op.Div(flux(op.internalConn));
df = Div(state.flux);
dms = Div(state_ms_it.flux);

clf; hold on
plot(df, '.');
plot(dms, 'o');
legend('Divergence (reference)', 'Divergence MS')

%% Use the class-based linear solver
% There is also a multiscale linear solver class for use with the ad-core
% framework, which uses the same internal interface.
mrstModule add ad-core
solver = MultiscaleVolumeSolverAD(CG, 'getSmoother', fn, ...
                                      'maxIterations', maxIterations,...
                                      'tolerance', tol, ...
                                      'useGMRES', true);

%% Solve a system
% We solve a linear system with the solveLinearSystem interface. Note that
% this generates a basis from A, which can then be reused for subsequent
% solves.
sol = solver.solveLinearSystem(A, b);
disp(solver.basis)

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
