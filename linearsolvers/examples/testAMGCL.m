%% Example demonstrating AMGCL on a few test problems
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
    layers = 1:5;
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
%% Compare with AGMG (if available)
try
    mrstModule add agmg
    tic();
    x_agmg = agmg(A, b, [], [], [], true);
    t_agmg = toc();
catch
    x_agmg = nan*b;
    t_agmg = nan;
end
%% Call AMGCL in AMG mode
tic();
x = callAMGCL(A, b, 'coarsening', 'ruge_stuben',...
                    'relaxation', 'spai0', ...
                    'preconditioner', 'amg', 'verbose', true);
t_amg = toc();
%% Call AMGCL in preconditioner mode
tic();
x = callAMGCL(A, b, 'relaxation', 'ilu0',...
                    'maxIterations', 1000, ...
                    'preconditioner', 'relaxation', 'verbose', true);
t_relax = toc();
%% Plot time taken
clf;
bar([t_amg, t_relax, t_agmg])
set(gca, 'XTicklabel', {'AMGCL-amg', 'AMGCL-relax' 'AGMG'});
ylabel('Solution time [s]');

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
