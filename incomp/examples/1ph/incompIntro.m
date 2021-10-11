%% Basic Flow-Solver Tutorial
% The purpose of this example is to give an overview of how to set up and
% use a standard two-point pressure solver to solve the single-phase
% pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\frac{K}{\mu}\nabla p,$$
%
% for a flow driven by Dirichlet and Neumann boundary conditions. Our
% geological model will be simple a Cartesian grid with anisotropic,
% homogeneous permeability.
%
% In this tutorial example, you will learn about:
%
% # the grid structure,
% # how to specify rock and fluid data,
% # the structure of the data-objects used to hold solution,
% # how to assemble and solve linear systems,
% # useful routines for visualizing and interacting with the grids and
% simulation results.

mrstModule add incomp

%% Define geometry
% Construct a Cartesian grid of size 10-by-10-by-4 cells, where each cell
% has dimension 1-by-1-by-1. Because our flow solvers are applicable for
% general unstructured grids, the Cartesian grid is here represented using
% <matlab:help('grid_structure') an unstructured format>, in which cells,
% faces, nodes, etc. are given explicitly.
nx = 10; ny = 10; nz = 4;
G = cartGrid([nx, ny, nz]);
display(G);

%%
% After the grid structure is generated, we plot the geometry.
clf, plotGrid(G); view(3)

%% Process geometry
% Having set up the basic structure, we continue to compute centroids and
% volumes of the cells and centroids, normals, and areas for the faces. For
% a Cartesian grid, this information can trivially be computed, but is
% given explicitly so that the flow solver is compatible with fully
% unstructured grids.
G = computeGeometry(G);


%% Set rock and fluid data
% The only parameters in the single-phase pressure equation are the
% permeability $K$ and the fluid viscosity $\mu$. We set the permeability
% to be homogeneous and anisotropic
%
% $$ K = \left(\begin{array}{ccc}
%      1000 & 0  & 0 \\ 0 & 100 & 0 \\ 0 & 0 & 10 \end{array}\right) $$
%
% The viscosity is specified by saying that the reservoir is filled with a
% <matlab:help('initSinglefluid') single fluid>, for which de default
% viscosity value equals unity. Our flow solver is written for a general
% incompressible flow and requires the evaluation of a total mobility,
% which is provided by the <matlab:help('fluid') fluid object>.
gravity off
rock = makeRock(G, [1000, 100, 10].* milli*darcy(), 1);
fluid     = initSingleFluid('mu' ,    1*centi*poise     , ...
                            'rho', 1014*kilogram/meter^3);

%% Initialize reservoir simulator
% To simplify communication among different flow and transport solvers, all
% unknowns are collected in a structure. Here this structure is initialized
% with uniform initial reservoir pressure equal 0 and (single-phase)
% saturation equal 0.0 (using the default behavior of
% <matlab:help('initResSol') initResSol>)
resSol = initResSol(G, 0.0);
display(resSol);

%% Impose Dirichlet boundary conditions
% Our flow solvers automatically assume no-flow conditions on all outer
% (and inner) boundaries; other type of boundary conditions need to be
% specified explicitly.
%
% Here, we impose Neumann conditions (flux of 1 m^3/day) on the global
% left-hand side. The fluxes must be given in units of m^3/s, and thus we
% need to divide by the number of seconds in a day (<matlab:help('day')
% day()>).  Similarly, we set Dirichlet boundary conditions p = 0 on the
% global right-hand side of the grid, respectively. For a single-phase
% flow, we need not specify the saturation at inflow boundaries. Similarly,
% fluid composition over outflow faces (here, right) is ignored by pside.
bc = fluxside([], G, 'LEFT',  1*meter^3/day());
bc = pside   (bc, G, 'RIGHT', 0);
display(bc);

%% Construct linear system
% Compute transmissibilities based on input grid and rock properties and
% use this to set up the linear system.
T = computeTrans(G, rock, 'Verbose', true);

% Solve linear system construced from T and bc to obtain solution for flow
% and pressure in the reservoir. The <matlab:help('incompTPFA')
% option> 'MatrixOutput=true' adds the system matrix A to resSol to enable
% inspection of the matrix.
resSol = incompTPFA(resSol, G, T, fluid, ...
                   'bc', bc, 'MatrixOutput', true);
display(resSol);

%% Inspect results
% The  |resSol| object contains the matrix used to solve the system.
spy(resSol.A);

%%
% We then plot convert the computed pressure to <matlab:help('units') unit>
% <matlab:help('barsa') 'bar'> before plotting result.
clf
plotCellData(G, convertTo(resSol.pressure(1:G.cells.num), barsa()), ...
             'EdgeColor', 'k');
title('Cell Pressure [bar]')
xlabel('x'), ylabel('y'), zlabel('Depth');
view(3); shading faceted; camproj perspective; axis tight;
colorbar

%% Copyright notice

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
