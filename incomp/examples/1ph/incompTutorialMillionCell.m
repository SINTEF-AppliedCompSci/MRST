%% Simulation of a Mega-cell Model
% The purpose of this example is to show how one can setup MRST so that it
% is capable of solving models with multi-million cells. To this end, we
% will use the two-point flux-approximation scheme combined with a highly
% efficient algebraic multigrid method (AGMG) to solve the single-phase
% pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\frac{K}{\mu}\nabla p,$$
%
% for a flow driven by Dirichlet and Neumann boundary conditions. Our
% geological model will be simple a Cartesian grid with anisotropic,
% homogeneous permeability.

mrstModule add incomp agmg linearsolvers

%% Setup iterative linear solver
% When working with large models, one cannot use the standard MLDIVIDE
% ('\') solver in MATLAB. Here, we use the AGMG or AMGCL algebraic
% multigrid solver.
if norm(callAMGCL(speye(3), [ 1 ; 2 ; 3 ]) - [ 1 ; 2 ; 3 ]) < 1.0e-8
    linsolver = @(A,b) callAMGCL(A,b);
    disp('Using AMGCL as linear solver');
elseif exist('agmg', 'file') && ...
      norm(agmg(speye(3), [ 1 ; 2 ; 3 ]) - [ 1 ; 2 ; 3 ]) < 1.0e-8
    linsolver = @(A,b) agmg(A,b,1);
    disp('Using AGMG as linear solver');
else
    error('This example requires the AGMG or AMGCL linear solvers');
end

%% Define geometry
% Construct a Cartesian grid of size 200-by-100-by-100 (=2,000,000) cells
% in a box of dimensions 10-by-10-by-4 metres. We have been able to run
% this problem with a peak memory use of slightly less than 4GB.
fprintf('Constructing grid ... ');
nx = 200; ny = 100; nz = 100;
G = cartGrid([nx, ny, nz], [10, 10, 4]);
fprintf('done\n');

%% Process geometry
% The computation of geometric information (centroids and volumes of the
% cells and centroids, normals, and areas for the faces) was accelerated
% significantly in MRST R2016a, but is still faster if one uses the
% C-accelerated routine instead of the standard call to computeGeometry(G)
mrstModule add libgeometry incomp
fprintf('Computing geometry ...');
G = mcomputeGeometry(G);
fprintf('done\n');

%% Set rock and fluid data
% The only parameters in the single-phase pressure equation are the
% permeability $K$ and the fluid viscosity $\mu$.
fprintf('Setting up rock and fluid parameters ...');
rock = makeRock(G, [1000, 100, 10].* milli*darcy(), 1);
T         = computeTrans(G, rock);

gravity off;
fluid     = initSingleFluid('mu' ,    1*centi*poise     , ...
                            'rho', 1014*kilogram/meter^3);
fprintf('done\n');

%% Set boundary conditions
% Here, we impose Neumann conditions (flux of 1 m^3/day) on the global
% left-hand side. Similarly, we set Dirichlet boundary conditions p = 0 on
% the global right-hand side of the grid. For a single-phase
% flow, we need not specify the saturation at inflow boundaries and the
% fluid composition over outflow faces (here, right) is ignored by pside.
fprintf('Initializing and setting boundary conditions ...');
resSol = initResSol(G, 0.0);
bc     = fluxside([], G, 'LEFT',  1*meter^3/day());
bc     = pside   (bc, G, 'RIGHT', 0);
fprintf('done\n');

%% Solve the linear system
% Solve linear system construced from T and bc to obtain solution for flow
% and pressure in the reservoir.
fprintf('Solving linear system ...');
resSol = incompTPFA(resSol, G, T, fluid, 'bc', bc, ...
                    'LinSolve', @(A,b) linsolver(A,b));
fprintf('done\n');

fprintf('Plotting results ...');
clf
plotCellData(G, convertTo(resSol.pressure(1:G.cells.num), barsa()),'EdgeColor','none');
title('Cell Pressure [bar]')
xlabel('x'), ylabel('y'), zlabel('Depth');
view(3); camproj perspective; axis tight;
colorbar
fprintf('done\n');

%% Copyright notice

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
