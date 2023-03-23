%% Quarter five spot
% The purpose of this example is to give an overview of how to set up and
% use the single-phase TPFA pressure solver to solve the pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\frac{K}{\mu}\nabla p,$$
%
% with no-flow boundary conditions and two source terms at diagonally
% opposite corners of a 2D Cartesian grid. This setup mimics a quarter
% five-spot well pattern, which is a standard test in reservoir simulation.
% In addition to computing pressure, we will also compute streamlines from
% the flux field as well as the corresponding time-of-flight, i.e., the
% time it takes for a neutral particle to travel from a fluid inlet (source
% term, well, inflow boundary, etc) to a given point in the reservoir.
mrstModule add incomp

%% Set up grid and petrophysical data
% We use a Cartesian grid of size nx-by-ny with homogeneous petrophysical
% data: permeability of 100 mD and porosity of 0.2.
[nx,ny] = deal(32);
G = cartGrid([nx,ny],[500,500]);
G = computeGeometry(G);
rock = makeRock(G, 100*milli*darcy, .2);

%% Compute half transmissibilities
% All we need to know to develop the spatial discretization is the reservoir
% geometry and the petrophysical properties. This means that we can compute
% the half transmissibilities without knowing any details about the fluid
% properties and the boundary conditions and/or sources/sinks that will
% drive the global flow:
hT = simpleComputeTrans(G, rock);


%% Fluid model
% When gravity forces are absent, the only fluid property we need in the
% incompressible, single-phase flow equation is the viscosity. However, the
% flow solver is written for general incompressible flow and requires the
% evaluation of a fluid object that can be expanded to represent more
% advanced fluid models. Here, however, we only use a simple fluid object
% that requires a viscosity and a density (the latter is needed when gravity
% is present)
gravity reset off
fluid = initSingleFluid('mu' , 1*centi*poise, ...
                        'rho', 1014*kilogram/meter^3);
display(fluid)

%% Add source terms
% To drive the flow, we will use a fluid source at the SW corner and a fluid
% sink at the NE corner. The time scale of the problem is defined by the
% strength of the source term. In our case, we set the source terms such
% that a unit time corresponds to the injection of one pore volume of
% fluids. All flow solvers in MRST automatically assume no-flow conditions
% on all outer (and inner) boundaries if no other conditions are specified
% explicitly.
pv  = sum(poreVolume(G,rock))/(10*year);
src = addSource([], 1, -pv);
src = addSource(src, G.cells.num, pv);
display(src)

%% Construct reservoir state object
% To simplify communication among different flow and transport solvers, all
% unknowns (reservoir states) are collected in a structure. Strictly
% speaking, this structure need not be initialized for an incompressible
% model in which none of the fluid properties depend on the reservoir
% states.  However, to avoid treatment of special cases, MRST requires that
% the structure is initialized and passed as argument to the pressure
% solver. We therefore initialize it with a dummy pressure value of zero and
% a unit fluid saturation since we only have a single fluid
state = initResSol(G, 0.0, 1.0);
display(state)

%% Solve pressure and show the result
% To solve for the pressure, we simply pass the reservoir state, grid model,
% half transsmisibilities, fluid model, and driving forces to the flow
% solver that assembles and solves the incompressible flow equation.
state = simpleIncompTPFA(state, G, hT, fluid, 'src', src);
display(state)

clf,
plotCellData(G, state.pressure);
plotGrid(G, src.cell, 'FaceColor', 'w');
axis equal tight; colormap(jet(128));

%% Trace streamlines
% To visualize the flow field, we show streamlines. To this end, we will
% use Pollock's method which is implemented in the 'streamlines' add-on
% module to MRST. Starting at the midpoint of all cells along the NW--SE
% diagonal in the grid, we trace streamlines forward and backward using
% 'pollock' and plot them using Matlab's builtin 'streamline' routine
mrstModule add streamlines;
seed = (nx:nx-1:nx*ny).';
Sf = pollock(G, state, seed, 'substeps', 1);
Sb = pollock(G, state, seed, 'substeps', 1, 'reverse', true);
hf=streamline(Sf);
hb=streamline(Sb);
set([hf; hb],'Color','k');

%% Compute time-of-flight
% Finally, we will compute time-of-flight, i.e., the time it takes a
% neutral particle to travel from the fluid source to a given point in the
% reservoir. Isocontours of the time-of-flight define natural time lines in
% the reservoir, and to emphasize this fact, we plot the time-of-flight
% using only a few colors.
mrstModule add diagnostics
tof = computeTimeOfFlight(state, G, rock, 'src', src);

clf,
plotCellData(G, tof);
plotGrid(G,src.cell,'FaceColor','w');
axis equal tight;
colormap(jet(16)); caxis([0,1]);

%% Visualize high-flow and stagnant regions
% We can also compute the backward time-of-flight, i.e., the time it takes
% a neutral particle to travel from a given point in the reservoir to the
% fluid sink. The sum of the forward and backward time-of-flights give the
% total travel time in the reservoir, which can be used to visualize
% high-flow and stagnant regions
tofb   = computeTimeOfFlight(state, G, rock, 'src', src, 'reverse', true);

clf,
plotCellData(G, tof+tofb);
plotGrid(G,src.cell,'FaceColor','w');
axis equal tight;
colormap(jet(128));

%% Exercises
% Compute a quarter five-spot setup using
% - perturbed grids, e.g., as computed by 'twister'
% - heterogeneous permeability and porosity, e.g., from a Karman-Cozeny
% relationship or as subsamples from SPE10

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
