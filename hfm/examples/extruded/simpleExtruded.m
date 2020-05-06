%% Single phase 2.5D problem
% The matrix and fractures are represented in 3D by extruding their
% corresponding 2D grids.

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add coarsegrid;      % functionality for coarse grids
mrstModule add ad-core;         % NNC support for coarse grids
mrstModule add msrsb;           % MsRSB solvers
mrstModule add mrst-gui;        % plotting routines
mrstModule add incomp;          % Incompressible fluid models
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path

%% Grid and fracture lines
% We start by constructing a 2D Cartesian grid comprising 50-by-50 cells,
% where each cell has dimension 1-by-1 m^2. Next, define 2 fracture lines
% in the form of a 'x' through the centre of the domain. We also define the
% number of layers in the extruded 3D grid and the indices of layers
% through which the fracture grid would extend. In this example, we extrude
% the matrix grid by 50 layers resulting in 50-by-50-by-50 fine cells. The
% fracture grid is extruded from layer 11 through 40.
celldim = [50 50];
G = cartGrid(celldim);
G = computeGeometry(G);
fl = [5 5 45 45; 5 45 45 5];
layers = 50;
flayers = 11:40;

%% Process fracture lines
% Using the input fracture lines, we identify independent fracture networks
% comprising of connected lines. In this example, there is only 1 fracture
% network consisting of 2 fracture lines. We also identify the fine-cells
% in the matrix containing these fractures. Fracture aperture is set to
% 0.04 meters. The matrix grid and fracture lines are plotted.
dispif(mrstVerbose, 'Processing user input...\n\n');
a = 1/25;
[G,fracture] = processFracture2D(G,fl); fracture.aperture = a;
figure;
plotFractureLines(G,fracture);
axis tight; box on

%% Compute CI and construct fracture grid
% For each matrix block containing a fracture, we compute a fracture-matrix
% conductivity index (CI) in accordance with the hierarchical fracture
% model (a.k.a. embedded discrete fracture model). Following that, we
% compute a fracture grid where the fracture cell size is defined to be
% 0.5 m. Next, the fracture grid is plotted on top of the matrix grid.
dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 0.5; cell_size = 0.5; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); box on

%% Make layered grid
% Using makeLayers, we extrude the matrix grid and each fracture grid along
% the z-direction. The fracture grid is extruded and readjusted in
% accordance with 'flayers'.
Gl = makeLayers(G,layers,flayers);

%% Set rock properties in fracture and matrix
% For this example, we will generate the matrix porosity as a Gaussian
% field. To get a crude approximation to the permeability-porosity
% relationship, we assume that our medium is made up of uniform spherical
% grains of diameter dp = 10 m, for which the specic surface area is Av =
% 6 = dp. With these assumptions, using the Carman Kozeny relation, we can
% then calculate the isotropic matrix permeability (K). The rock properties
% are then plotted. Fracture permeability is set to 1000 Darcy with 50%
% porosity in each fracture grid cell. The rock properties are then
% plotted.
p = gaussianField(Gl.cartDims, [0.2 0.4], [9 5 5], 3.5);
K = p.^3.*(1e-5)^2./(0.81*72*(1-p).^2);

Gl.rock = makeRock(Gl, K(:), p(:));
K_frac = 10000; % Scaling factor = K_frac/K_mat D
poro_frac = 0.5;
Gl = makeRockFrac(Gl, K_frac, 'permtype','homogeneous','porosity',poro_frac);

clf; plotToolbar(Gl,Gl.rock);
colormap(jet); view (-135,30);
axis tight equal off

%% Define fluid properties
% Define a single fluid of viscosity 1 cP and density 1000 kg/m3.
fluid = initSingleFluid('mu' , 1*centi*poise, ...
    'rho', 1000*kilogram/meter^3);

%% Global grid with NNC's and corresponding transmissibility
% In this section, we use the function makeNNCextruded to combine the
% layered fracture and matrix grid structures into a single grid structure.
% Additionally, we assign a 'non-neighbouring connection (NNC)' status to
% every fracture-matrix connection. For example, if fracture cell 'i' is
% located in the matrix cell 'j', then the connection i-j is defined as an
% NNC. To compute the flux between these elements, we compute a
% transmissibility for each NNC using the CI's computed earlier. Vector T
% contains the transmissibility for each face in the combined grid and each
% NNC.
[Gl,T] = makeNNCextruded(G,Gl,F,fracture,flayers);
G = Gl; clear Gl

%% Add Wells
% We will include two horizontal wells, a rate-controlled injector well
% through the bottom left corner of the grid and a producer, controlled by
% bottom-hole pressure, through the top right section of the grid. Wells
% are described using a Peaceman model, giving an extra set of equations
% that need to be assembled, see simpleWellExample.m for more details on
% the Peaceman well model.
[nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
cellinj = nx*ny*(nz-1)+1:nx:nx*ny*nz;
cellprod = nx:nx:nx*ny;

W = [];
W   = addWell(W, G.Matrix, G.Matrix.rock, cellinj, 'Type', 'rate',...
    'Val', 1000/day, 'Sign',1, 'Comp_i', 1, 'Name', 'Injector');
W   = addWell(W, G.Matrix, G.Matrix.rock, cellprod, ...
    'Type', 'bhp', 'Val', 50*barsa(), 'Sign', -1, 'Comp_i', 1, 'Name', 'Producer');

%% Initialize state variables
% Once the transmissibilities are computed, we can generate the
% transmissiblity matrix 'A' using the 'two-point flux approximation'
% scheme and initialize the solution structure.
dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 0);
state.wellSol = initWellSol(W, 0);
[A,q] = getSystemIncompTPFA(state, G, T, fluid, 'use_trans', true);

%% Solve fine-scale problem
% The fine-scale pressure solution is computed using the boundary
% conditions provided and the transmissiblity matrix computed earlier.
dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state_fs = incompTPFA(state, G, T, fluid,  ...
    'Wells', W, 'MatrixOutput', true, 'use_trans',true);

%% Setup multiscale grids
% Next, we define a 10-by-10-by-10 matrix coarse grid such that each coarse
% block contains 5-by-5-by-5 fine cells. Each of the two fracture planes is
% partitioned into 3-by-3 coarse blocks. In total, the fracture grid has 18
% degrees of freedom at coarse scale giving a coarsening ratio of ~333.
% Additionally, we also define the support regions for the fracture and
% matrix basis functions. Fracture support region is defined based on a
% topological distance based algorithm. The matrix and fracture coarse
% grids are plotted in the next section.
G.type{1,1} = 'layered';

% Partition matrix
coarseDims = [10 10 10];
pm = partitionMatrix(G, 'coarseDims', coarseDims, 'use_metis', false);
CGm = getRsbGridsMatrix(G, pm, 'fullyCoupled', false, 'Wells', W);

% Partition fracture
nw = fracture.network;
coarseDimsF = [3 3];
p  = partitionFracture(G, pm, nw, 'partition_frac'   , true   , ...
    'use_metisF'       , false  , ...
    'coarseDimsF'      , coarseDimsF );

p = processPartition(G,compressPartition(p));
pf = p(G.Matrix.cells.num+1:end)-max(p(1:G.Matrix.cells.num));

% Coarse Grids
CG = generateCoarseGrid(G, p);

% Add centroids / geometry information on coarse grid
CG = coarsenGeometry(CG);
Gf = assembleFracGrid(G);
CGf = generateCoarseGrid(Gf, pf);
CGf = coarsenGeometry(CGf);

% Support Regions
[CG,CGf] = storeFractureInteractionRegion(CG, CGf, CGm, ...
    'excludeBoundary' , false , ...
    'removeCenters'   , false , ...
    'fullyCoupled'    , false );


%% Compute basis functions
% Using the transmissibility matrix 'A' we compute the basis functions for
% each fracture and matrix coarse block using the restriction smoothed
% basis method. Note that the matrix 'A' does not contain any source terms
% or boundary conditions. They are added to the coarse linear system when
% computing the multiscale pressure in the next section.
dispif(mrstVerbose, 'Computing basis functions...\n\n');
basis_sb = getMultiscaleBasis(CG, A, 'type', 'msrsb');
clf; plotToolbar(G,basis_sb.B,'filterzero',true); view(-135,30)
plotGrid(CG, 'FaceColor', 'none', 'linewidth', 1);
axis tight; colormap(jet); colorbar;
title('Basis Functions in the matrix');

%% Compute multiscale solution
state_ms = incompMultiscale(state, CG, T, fluid, basis_sb, 'Wells', W, ...
    'use_trans',true);

 %% Plot results
figure;
plotToolbar(G, state_fs.pressure)
colormap jet; colorbar
view(-135,30)
axis tight off
title('Fine scale')

figure;
plotToolbar(G, state_ms.pressure)
colormap jet; colorbar
view(-135,30)
axis tight off
title('F-MsRSB')

L1 = abs(state_ms.pressure-state_fs.pressure)./state_fs.pressure;
figure;
plotToolbar(G, L1)
colormap jet; colorbar
view(-135,30)
axis tight off
L1_eq = '$$| P_i^{fs}-P_i^{ms} | / P_i^{fs} $$';
title(L1_eq,'interpreter','latex');

%% Solve using MS-ILU and MS-GMRES
% As seen in the plots above, the multiscale method is able to capture the
% general trend of the solution, but can have large pointwise errors. To
% improve the approximation, we can compute iterative multiscale solutions
% using ILU and GMRES preconditioners.
fn = getSmootherFunction('type', 'ilu');

[~,report] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'Wells', W, 'use_trans', true, 'tolerance', 1e-8, 'iterations', 100,...
    'useGMRES', false, 'reconstruct', true, 'getSmoother', fn);

[~,report2] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'Wells', W, 'use_trans', true, 'tolerance', 1e-8, 'iterations', 100,...
    'useGMRES', true, 'reconstruct', true, 'getSmoother', fn);


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