%% Pressure Solver: Simple Corner-Point Grid with Linear Pressure Drop
% Herein we will solve the single-phase pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\frac{K}{\mu}\nabla p,$$
%
% for a simple corner-point grid model with isotropic, lognormal
% permeability.
%
% The purpose of this example is to demonstrate the two-point pressure
% solver applied to a case with corner-point grids given as from an input
% stream in the industry-standard Eclipse(TM) format.
%

mrstModule add incomp

%% Generate the corner-point grid
% The corner-point grid is generated using a standard MATLAB(R) meshgrid
% which is then transformed to make sloping pillars and wavy layers. The
% corner-point grid is represented as a GRDECL structure, which is the same
% structure as is returned from readGRDECL.  If an extra parameter is
% passed to the simpleGrdecl function, it adds a fault to the model. 
nc = [20, 20, 5];
grdecl = simpleGrdecl(nc);        % No fault in model

%%
% Then we process the GRDECL structure and build up an unstructured grid.
G = processGRDECL(grdecl); clear grdecl

%%
% After the grid structure is generated, we plot the geometry.
clf, plotGrid(G);
title('Corner-point grid')
view(3), camproj perspective, axis tight, camlight headlight

%% Compute geometry information
% Having set up the basic structure, we continue to compute centroids and
% volumes of the cells and centroids, normals, and areas for the faces.
G = computeGeometry(G, 'Verbose', true);


%% Make rock and fluid data
% We generate a lognormal and layered permeability and specify that we are
% working with a single-phase fluid; type <matlab:helpt('initSingleFluid')
% "help initSingleFluid"> to see default values for density and viscosity.
% Our flow solver is constructed for a general incompressible flow and thus
% requires the evaluation of a total mobility function, which in our case
% equals one in the whole domain.
perm = convertFrom(logNormLayers(nc, [100, 400, 50]), milli*darcy());
rock = makeRock(G, perm, 1);
fluid     = initSingleFluid('mu' ,    1*centi*poise     , ...
                            'rho', 1014*kilogram/meter^3);

%%
% Plot the logarithm of the layered permeability.
cla,
plotCellData(G, log10(rock.perm(:)));
title('Log_{10} of permeability')
camproj perspective, axis tight off, camlight headlight

%% Add boundary conditions
% Pressure is set of one bar on the west-side boundary and zero on the
% east-side boundary.  No-flow conditions are imposed automatically at all
% other (outer) boundaries.
westFaces = find(G.faces.centroids(:,1) == 0);
bc = addBC([], westFaces, 'pressure',        ...
           repmat(1*barsa(), [numel(westFaces), 1]));

xMax = max(G.faces.centroids(:,1));
eastFaces = find(G.faces.centroids(:,1) == xMax);
bc = addBC(bc, eastFaces, 'pressure',        ...
           zeros(numel(eastFaces), 1));

%%
% Then we plot the grid, coloring the faces on which we have imposed
% boundary conditions.
cla,
plotGrid(G, 'FaceColor', 'none');
plotFaces(G, westFaces, 'r');
plotFaces(G, eastFaces, 'b');
title('Boundary conditions')
camproj perspective, axis tight off, camlight headlight


%%  Assemble and solve system
% Finally, we compute transmissibilities, construct two-point
% discretization and solve the corresponding linear equations.
rSol = initResSol(G, 0);
T    = computeTrans(G, rock, 'Verbose', true);
rSol = incompTPFA(rSol, G, T, fluid, 'MatrixOutput', true, 'bc', bc);

%%
% We plot the results: cell pressure is converted to unit 'bar' and fluxes
% to unit m^3/day when plotting results.
subplot('Position',[0.01 0.25 0.48 0.5]),
   plotCellData(G, convertTo(rSol.pressure(1:G.cells.num), barsa()));
   title('Pressure [bar]')
   view(3), camproj perspective, axis tight off, camlight headlight

subplot('Position',[0.51 0.25 0.48 0.5]),
   cellNo  = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
   plotCellData(G, sqrt(accumarray(cellNo, ...
      abs(convertTo(faceFlux2cellFlux(G, rSol.flux), meter^3/day)))));
   title('Sqrt - Flux intensity [m^3/day]')
   view(3), camproj perspective, axis tight off, camlight headlight

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

