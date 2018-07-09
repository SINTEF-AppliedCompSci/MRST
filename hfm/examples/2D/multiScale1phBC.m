%% Introduction to the F-MsRSB Solver
% In this example, we will introduce you to the multiscale restriction
% smoothed basis (MsRSB) method for computing flow in embedded fracture
% models. To this end, we consider a 2D single-phase example with two
% intersecting fractures and Dirichelet boundary conditions. The flow
% problem is solved both by a fine-scale and a multiscale solver.
%
% Notice that you need to have Metis installed to get this example to work.
% To get Metis working, you also need to set the global variable METISPATH.
% This can be done in your 'startup_user.m' file.

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add coarsegrid;      % functionality for coarse grids
mrstModule add ad-core;         % NNC support for coarse grids
mrstModule add msrsb;           % MsRSB solversop
mrstModule add mrst-gui;        % plotting routines
mrstModule add incomp;          % Incompressible fluid models
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path
    
%% Grid and fracture lines
% Construct a Cartesian grid comprising 90-by-90 cells, where each cell has
% dimension 0.1-by-0.1 m^2. Define 2 fracture lines, 5 m in length, in
% the form of a '+' in the centre of the domain.

celldim = [90 90];
physdim = [9 9];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

fl = [20, 45, 70, 45; 
      45, 20, 45, 70]./10 ; % fractures lines in [x1 y1 x2 y2] format.

%% Process fracture lines
% Using the input fracture lines, we identify independent fracture networks
% comprising of connected lines. In this example, there is only 1 fracture
% network consisting of 2 fracture lines. We also identify the fine-cells
% in the matrix containing these fractures. Fracture aperture is set to
% 0.04 meters. The matrix grid and fracture lines are plotted.

dispif(mrstVerbose, 'Processing user input...\n\n');
[G,fracture] = processFracture2D(G,fl);
fracture.aperture = 1/25; % Fracture aperture
figure;
plotFractureLines(G, fracture);
axis tight; 
box on

%% Compute CI and construct fracture grid
% For each matrix block containing a fracture, we compute a fracture-matrix
% conductivity index (CI) in accordance with the hierarchical fracture
% model (a.k.a. embedded discrete fracture model). Following that, we
% compute a fracture grid where the fracture cell size is defined to be
% 0.05 m. Next, the fracture grid is plotted on top of the matrix grid.

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 0.05; cell_size = 0.05; % minimum and average cell size.
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

%% Initialize state variables
% Once the transmissibilities are computed, we can generate the
% transmissiblity matrix 'A' using the 'two-point flux approximation'
% scheme and initialize the solution structure.

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 5*barsa);

% Get A matrix without source
[A,q] = getSystemIncompTPFA(state, G, T, fluid, 'use_trans', true);

%% Setup multiscale grids
% Next, we define a 9-by-9 matrix coarse grid and 4 coarse blocks (or
% coarse degrees of freedom) in the fracture. Additionally, we also define
% the support regions for the fracture and matrix basis functions. The
% matrix and fracture coarse grids are plotted.

dispif(mrstVerbose, 'Defining coarse grids and interaction regions...\n\n');

coarseDims = [9 9]; 
dof_frac = 4; % Number of coarse blocks in the fracture grid
[CG, CGf] = getRsbGridsHFM(G, fracture.network, 'coarseDims', coarseDims,...
            'dof_frac',dof_frac);

clf; plotFractureCoarseGrid2D(G,CG.partition,F)

%% Add BC
% Set boundary condition P = 10 bars on the left face and P = 1 bar on the
% right face of the domain.

bc = [];
bc  = pside(bc, G, 'LEFT', 10*barsa);
bc  = pside(bc, G, 'RIGHT', 1*barsa);

%% Incompressible fine-scale solver
% The fine scale pressure solution is computed using the boundary
% conditions provided and the transmissiblity matrix computed earlier.

dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state_fs = incompTPFA(state, G, T, fluid,  ...
            'bc', bc, 'MatrixOutput', true, 'use_trans',true);

%% Compute basis functions
% Using the transmissibility matrix 'A' we compute the basis functions for
% each fracture and matrix coarse block using the restriction smoothed
% basis method. Note that the matrix 'A' does not contain any source terms
% or boundary conditions. They are added to the coarse linear system when
% computing the multiscale pressure in the next section.

dispif(mrstVerbose, 'Computing basis functions...\n\n');
basis_sb = getMultiscaleBasis(CG, A, 'type', 'rsb');
clf; plotToolbar(G,basis_sb.B,'filterzero',true);
plotGrid(CG,'FaceColor','none');
axis tight; colorbar; 
title('Basis Functions in the matrix');

%% Compute multiscale solution

dispif(mrstVerbose, 'Computing multiscale solution...\n\n');
[state_ms,~] = incompMultiscale(state, CG, T, fluid, basis_sb,...
               'bc', bc,'use_trans',true);

%% Solve using MS-ILU and MS-GMRES
% Compute an iterative multiscale solution using ILU and GMRES
% preconditioners.

fn = getSmootherFunction('type', 'ilu');

[~,report] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'bc', bc, 'use_trans',true, 'tolerance', 1e-8, 'iterations', 100,...
    'useGMRES', false, 'reconstruct', true, 'getSmoother', fn);

[~,report2] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'bc', bc, 'use_trans',true, 'tolerance', 1e-8, 'iterations', 100,...
    'useGMRES', true, 'reconstruct', true, 'getSmoother', fn);

%% Plot results and convergence
figure; colormap jet(25)
plotCellData(G, state_fs.pressure,'EdgeColor','none')
line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
view(90, 90); colorbar, cx=caxis();
axis tight off
title('Fine scale')

figure; colormap jet(25)
plotCellData(G, state_ms.pressure,'EdgeColor','none')
line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
view(90, 90); colorbar, caxis(cx);
axis tight off
title('F-MsRSB')

figure; colormap jet
L1 = abs(state_ms.pressure-state_fs.pressure)./state_fs.pressure;
plotCellData(G, L1,'EdgeColor','none')
line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
view(90, 90); colorbar
axis tight off
L1_eq = '$$ \frac{| P_i^{fs}-P_i^{f-msrsb} | }{ P_i^{fs}} $$';
title(L1_eq,'interpreter','latex');

%% Plot convergence
figure;
semilogy(report.resvec, ':+'); hold on;
semilogy(report2.resvec, 's-');
legend('Convergence of ILU(0)', 'Convergence of GMRES');
title('Iterative convergence');

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