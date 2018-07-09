%% Single-Phase Problem Demonstrating the Hierarchical Fracture Model
% In this first introductory example to the HFM module, we consider a
% single-phase 2D example with Dirichlet boundary conditions and a
% horizontal central fracture in the centre. The example shows the effect
% of a single fracture on the reservoir pressure distribution. Moreover,
% the pressure solution obtained using an embedded discrete or hierarchical
% fracture model is compared to the results obtained from a fully resolved
% simulation where the fracture and matrix grid blocks are of the same
% size.

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add incomp;          % Incompressible fluid models
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path

%% Grid and fracture lines
% Construct a Cartesian grid comprising 90-by-90 cells, where each cell has
% dimension 0.1-by-0.1 m^2. Define a horizontal fracture line, 5 m in
% length, in the centre of the domain by supplying its end points (2, 4.5)
% and (7, 4.5) as a 1-by-4 array of the format [x1 y1 x2 y2].

celldim = [90 90];
physdim = [9 9];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

fl = [20, 45, 70, 45]./10 ; % fractures lines in [x1 y1 x2 y2] format.

%% Process fracture lines
% Using the input fracture lines, we identify independent fracture networks
% comprising of connected lines. In this example, there is only 1 fracture
% network consisting of 1 fracture line. We also identify the fine-cells in
% the matrix containing these fractures. Fracture aperture is set to 0.04
% meters. The matrix grid and fracture lines are plotted.

dispif(mrstVerbose, 'Processing user input...\n\n');
[G,fracture] = processFracture2D(G,fl);
fracture.aperture = 1/25; % Fracture aperture
figure;
plotFractureLines(G,fracture);
axis tight; 
box on

%% Compute CI and construct fracture grid
% For each matrix block containing a fracture, we compute a fracture-matrix
% conductivity index (CI) in accordance with the hierarchical fracture
% model (a.k.a. embedded discrete fracture model). Following that, we
% compute a fracture grid where the fracture cell size is defined to be
% 0.1 m. Next, the fracture grid is plotted on top of the matrix grid.

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 0.05; cell_size = 0.1; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); box on

%% Set rock properties in fracture and matrix
% Set the permeability (K) as 1 Darcy in the matrix and 10000 Darcy in the
% fractures.

dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
G.rock.perm = ones(G.cells.num,1) * darcy;
K_frac = 10000; % Darcy
G = makeRockFrac(G, K_frac);

%% Define fluid properties
% Define a single fluid of viscosity 1 cP and density 1000 kg/m3.

fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1000*kilogram/meter^3);

%% Define fracture connections as NNC and compute the transmissibilities
% In this section, we use the function defineNNCandTrans to combine the
% fracture and matrix grid structures into a single grid structure. In
% addition to that, we assign a 'non-neighbouring connection (NNC)' status
% to every fracture-matrix connection. For example, if fracture cell 'i' is
% located in the matrix cell 'j', then the connection i-j is defined as an
% NNC. To compute the flux between these elements, we compute a
% transmissibility for each NNC using the CI's computed earlier. Vector T
% contains the transmissibility for each face in the combined grid and each
% NNC.

[G,T] = defineNNCandTrans(G,F,fracture);

%% Define completely resolved grid
% We now define the fully resolved grid where the matrix and fracture cells
% are of the same size. Rock properties for this new grid are defined such
% that a 5 m long horizontal line of cells in the centre of the domain has
% a permeability of 10000 Darcy whereas the rest of the grid has a
% permeability of 1 Darcy. 

celldim = [225 225];
physdim = [9 9];
Gr = cartGrid(celldim, physdim);
Gr = computeGeometry(Gr);
Gr.rock.perm = ones(Gr.cells.num,1) * darcy;
fractureCells = 225*112+50:225*112+175;
Gr.rock.perm(fractureCells) = 10000*darcy;
clf; plotCellData(Gr,Gr.rock.perm,'EdgeColor', 'none'); axis equal tight

%% Initialize state variables
% Here, we initialize the solution structure for the two grids.

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 5*barsa);
state_r  = initResSol (Gr, 5*barsa);

%% Add BC
% Set boundary condition P = 10 bars on the left face and P = 1 bar on the
% right face of the domain.

bc = [];
bc  = pside(bc, G, 'LEFT', 10*barsa);
bc  = pside(bc, G, 'RIGHT', 1*barsa);

bcr = [];
bcr  = pside(bcr, Gr, 'LEFT', 10*barsa);
bcr  = pside(bcr, Gr, 'RIGHT', 1*barsa);

%% Compute pressure
% The fine scale pressure solution is computed for both the grids using the
% boundary conditions provided and the transmissiblity matrix computed
% earlier.

dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state = incompTPFA(state, G, T, fluid, 'bc', bc, 'use_trans',true);
state_r = incompTPFA(state_r, Gr, computeTrans(Gr,Gr.rock), fluid, 'bc', bcr);

%% Plot results
% Use a reduced number of color levels to better compare and contrast the
% two solutions
figure
plotCellData(Gr, state_r.pressure,'EdgeColor','none')
line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
colormap jet(25)
view(0, 90); colorbar; axis equal tight, cx = caxis();
title('Fully resolved solution')

figure
plotCellData(G, state.pressure,'EdgeColor','none')
line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
colormap jet(25)
view(0, 90); colorbar; axis equal tight, caxis(cx);
title('Hierarchical or Embedded discrete fracture model');

%% Add a contour plot
figure
contour(reshape(Gr.cells.centroids(:,1),Gr.cartDims),...
   reshape(Gr.cells.centroids(:,2),Gr.cartDims),...
   reshape(state_r.pressure,Gr.cartDims),25);
hold on
nc = prod(G.cartDims);
contour(reshape(G.cells.centroids(1:nc,1),G.cartDims),...
   reshape(G.cells.centroids(1:nc,2),G.cartDims),...
   reshape(state.pressure(1:nc),G.cartDims),25,'--','LineWidth',2);
line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
hold off
legend('Fully resolved','Embedded fractures');

% <html>
% <p><font size="-1">
% Copyright 2009-2016 TU Delft and SINTEF ICT, Applied Mathematics.
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
