%% Two-Phase Problem with two Intersecting Fractures
% Incompressible two-phase problem with water injection from the left
% boundary of a square domain and production at the opposite face. The grid
% contains 2 intersecting fractures. The flow problem is solved both by a
% fine-scale and a multiscale solver.
% 
% Notice that you need to have Metis installed to get this example to work.
% To get Metis working, you also need to set the global variable METISPATH.
% This can be done in your 'startup_user.m' file.


% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add coarsegrid;      % functionality for coarse grids
mrstModule add ad-core;         % NNC support for coarse grids
mrstModule add new-multiscale;  % MsRSB solvers
mrstModule add mrst-gui;        % plotting routines
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path
checkMATLABversionHFM;

%% Grid and fracture lines
% Construct a Cartesian grid comprising 99-by-99 cells, where each cell has
% dimension 1-by-1 m^2. Define 2 fracture lines in the form of a 'x'
% through the centre of the domain.

celldim = [99 99];
physdim = [99 99];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

fl = [ 22    22    77    77; 22    77    77    22];

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
plotFractureLines(G,fracture);
box on

%% Compute CI and construct fracture grid
% For each matrix block containing a fracture, we compute a fracture-matrix
% conductivity index (CI) in accordance with the hierarchical fracture
% model (a.k.a. embedded discrete fracture model). Following that, we
% compute a fracture grid where the fracture cell size is defined to be
% 0.5 m. Next, the fracture grid is plotted on top of the matrix grid.

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 0.1; cell_size = 0.5; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); box on

%% Set rock properties in fracture and matrix
% Set the permeability (K) as 1 Darcy in the matrix and 10000 Darcy in the
% fractures. Additionally, set the porosity of the matrix and fractures to
% 20% and 50%, respectively.

dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
G.rock.perm = ones(G.cells.num,1)*darcy;
G.rock.poro = 0.2*ones(G.cells.num, 1);
K_frac = 10000; % Darcy
poro_frac = 0.5;
G = makeRockFrac(G, K_frac, 'permtype','homogeneous','porosity',poro_frac);

%% Define fluid properties
% Define a two-phase fluid model without capillarity. The fluid model has
% the values:
%
% * densities: [rho_w, rho_o] = [1000 700] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 2] cP.
% * corey-coefficient: [2, 2] = [2 2].

fluid = initSimpleFluid('mu' , [   1,  2] .* centi*poise     , ...
    'rho', [1000, 700] .* kilogram/meter^3, ...
    'n'  , [   2,   2]);

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
xf = G.faces.centroids(:, 1);
left = find(abs(xf - min(xf)) < 1e-4);
right = find(abs(xf - max(xf)) < 1e-4);

bc = addBC(bc, left, 'pressure', 10*barsa, 'sat', [1 0]);
bc = addBC(bc, right, 'pressure', 1*barsa, 'sat', [0 1]);

%% Compute initial pressure

dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state_fs = incompTPFA(state, G, T, fluid,  ...
    'bc',bc, 'MatrixOutput', true, 'use_trans',true);

%% Compute basis functions
% Using the transmissibility matrix 'A' we compute the basis functions for
% each fracture and matrix coarse block using the restriction smoothed
% basis method. Note that the matrix 'A' does not contain any source terms
% or boundary conditions. They are added to the coarse linear system when
% computing the multiscale pressure in the next section.

dispif(mrstVerbose, 'Computing basis functions...\n\n');
basis_sb = getMultiscaleBasis(CG, A, 'type', 'rsb');
clf; plotToolbar(G,basis_sb.B);
axis tight; c = colormap(jet);
colormap(c); colorbar;
title('Basis functions plotted in the matrix');

%% Compute multiscale solution

dispif(mrstVerbose, 'Computing multiscale solution...\n\n');
[state_ms,~] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'bc', bc,'use_trans',true);

%% Plot initial pressure

figure;
plotToolbar(G, state_fs.pressure)
colormap jet
view(90, 90)
axis tight off
title('Initial Pressure: Fine scale')

figure;
plotToolbar(G, state_ms.pressure)
colormap jet
view(90, 90)
axis tight off
title('Initial Pressure: F-MsRSB')

%% Incompressible Two-Phase Flow
% We solve the two-phase system using a sequential splitting in which the
% pressure and fluxes are computed by solving the flow equation and then
% held fixed as the saturation is advanced according to the transport
% equation. Note that the flow equations are solved in the fine scale as
% well as on the coarse scale using an algebraic multiscale strategy. The
% multiscale flux field obtained at fine scale resolution is reconstructed
% to be conservative before solving the transport equation. This procedure
% is repeated for a given number of time steps (here we use 30 equally
% spaced time steps). The error introduced by this splitting of flow and
% transport can be reduced by iterating each time step until e.g., the
% residual is below a certain user-prescribed threshold (this is not done
% herein).

pv     = poreVolume(G,G.rock);
nt     = 30;
t200    = 2*(sum(pv)/sum(state_fs.flux(left)));
Time   = t200;
dT     = Time/nt;
dTplot = Time/3;
N      = fix(Time/dTplot);

state_ms.rhs = q;
t  = 0; plotNo = 1; hfs = 'Reference Saturation: '; hms = 'F-MsRSB Saturation: ';
dispif(mrstVerbose, 'Solving Transport...\n\n');
figure
B = basis_sb.B;
R = controlVolumeRestriction(CG.partition);
while t < Time,
    state_fs = implicitTransport(state_fs, G, dT, G.rock, fluid, 'bc', bc, 'Trans', T, 'verbose', true);
    state_ms = implicitTransport(state_ms, G, dT, G.rock, fluid, 'bc', bc, 'Trans', T);
    % Check for inconsistent saturations
    s = [state_fs.s(:,1); state_ms.s(:,1)];
    assert(max(s) < 1+eps && min(s) > -eps);
    
    % Update solution of pressure equation.
    state_fs  = incompTPFA(state_fs, G, T, fluid, 'bc', bc, 'use_trans', true);
    
    %-------------------------------Multiscale----------------------------%
    A = incompTPFA(state_ms, G, T, fluid, 'MatrixOutput', true, ...
        'use_trans',true); A = A.A;
    B = iteratedJacobiBasis(A, CG, 'interpolator', B);
    basis_sb = struct('B', B, 'R', R);
    state_ms = incompMultiscale(state_ms, CG, T, fluid, basis_sb,...
        'bc', bc, 'use_trans', true, 'reconstruct', true);
    %---------------------------------------------------------------------%
    
    % Increase time and continue if we do not want to plot saturations
    t = t + dT;
    if ( t < plotNo*dTplot && t < Time), continue, end
    
    % Plot saturation
    heading = [num2str(convertTo(t,day)),  ' days'];
    
    r = 0.01;
    subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.48]), cla
    plotCellData(G, state_fs.s(:,1), 'edgealpha', 0.1);
    colormap(flipud(jet))
    view(0,90), axis equal off, title([hfs heading])
    
    subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.48]), cla
    plotCellData(G, state_ms.s(:,1), 'edgealpha', 0.1);
    colormap(flipud(jet))
    view(0,90), axis equal off, title([hms heading])
    
    plotNo = plotNo+1;
    
end
