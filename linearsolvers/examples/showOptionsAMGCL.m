%% Example demonstrating different options for AMGCL
% Note that these choices mirror the names given in the AMGCL
% documentation: https://amgcl.readthedocs.io/en/latest/runtime.html
mrstModule add msrsb incomp coarsegrid spe10 linearsolvers
%% Setup problem
simple = false;
if simple
    % Create n by n grid
    n = 10;
    G = cartGrid([n, n, 1]);
    G = computeGeometry(G);
    rock = makeRock(G, 1, 1);
    W = [];
    W = addWell(W, G, rock, 1, 'type', 'rate', 'val', 1);
    W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 0);
else
    % Take part of SPE10 benchmark
    layers = 1;
    [G, W, rock] = getSPE10setup(layers);
end
% Set up linear system by calling the incompTPFA solver without a linear
% solver
fluid = initSingleFluid('rho', 0, 'mu', 1);
state = initResSol(G, 0);
state = incompTPFA(state, G, computeTrans(G, rock), fluid, 'Wells', W, 'matrixOutput', true, 'LinSolve', @(A, b) 0*b);
A = state.A;
b = state.rhs;
% Reduce to only wells
[A, b] = eliminateWellEquations(A, b, G.cells.num);

makeChoice = @(text, choices)  listdlg('PromptString',  text,...
                                       'SelectionMode', 'single',...
                                       'ListString',    choices);

%% Select preconditioner
% There are two options here. If relaxation is chosen, AMG will be used as
% a preconditioner. If relaxation is chosen, no multigrid will be used and
% the solver will be a simpler Krylov solver with relaxation as a
% preconditioner. The coarsening option chosen later in this example will
% not be used in this case.
preconditioners = {'amg', ...
                   'relaxation'};
index = makeChoice('Select solver', preconditioners);
precond = preconditioners{index};

%% Select solver
% Choices for Krylov acceleration
solvers = {'bicgstab', ... % Biconjugate gradient stabilized method
           'cg', ... % Conjugate gradient
           'bicgstabl', ... % Bicgstab variant
           'gmres', ... % Generalized minimal residual method
           'lgmres', ... % GMRES variant
           'fmgres', ... % Flexible-GMRES 
           'idrs'}; % Induced Dimension Reduction method IDR(s)
index = makeChoice('Select Krylov-solver', solvers);
solver = solvers{index};
%% Select coarsening and interpolation scheme for AMG
coarsening = {'ruge_stuben', ... % Ruge-stuben coarsening ("standard" coarsening)
              'smoothed_aggregation', ... % Aggregation with smoothed iterates
              'aggregation', ... % Aggregation - very fast coarsening due to very simple constant interpolation
              'smoothed_aggr_emin'}; % Aggregation with energy-minimizing smoothing
index = makeChoice('Select coarsening', coarsening);
coarsen = coarsening{index};
%% Select relaxation ("smoother")
relaxations = {'spai0', ... % Sparse approximate inverse, degree 0
               'spai1', ... % Sparse approximate inverse
               'gauss_seidel', ... % Gauss-seidel
               'ilu0', ... % Incomplete LU-factorization without fill-in
               'iluk', ... % Level-based ILU
               'ilut', ... % Thresholded ILU
               'damped_jacobi', ... % Damped Jacobi
               'chebyshev'}; % Chebyshev smoothing
index = makeChoice('Select relaxation', relaxations);
relax = relaxations{index};
%% Solve the system
tic();
[x, err, iter] = callAMGCL(A, b, ...
                    'coarsening', coarsen,...
                    'relaxation', relax, ...
                    'solver', solver, ...
                    'maxIterations', 100, ...
                    'preconditioner', precond, ...
                    'verbose', true);
t_amg = toc();
fprintf('%d iterations done in %f seconds. Final residual: %e\n', iter, t_amg, err);
%% Additional options
% There are many additional options which can be set through the AMGCL
% struct. Generally, these can be passed as keyword arguments to callAMGCL
% and other routines which rely on AMGCL.
% 
% Especially interesting options are
%  - ncycle (number of multigrid cycles)
%  - npre (number of presmootings)
%  - npost (number of postsmoothings)
% As well as options relating to the coarsening hierarchy:
%  - rs_eps_strong (strong threshold, ruge_stuben)
%  - aggr_eps_strong (strong threshold, aggregation)
str = getAMGCLMexStruct();

disp(str)

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
