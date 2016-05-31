%{
Single-phase 2D example with Dirichlet boundary conditions demonstrating
the use of the HFM module. The flow problem is solved both by a fine-scale
and a multiscale solver.

Notice that you need to have Metis installed to get this example to work.
To get Metis working, you also need to set the global variable METISPATH.
This can be done in your 'startup_user.m' file.
%}

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add coarsegrid;      % functionality for coarse grids
mrstModule add new-multiscale;  % MsRSB solvers
mrstModule add mrst-gui;        % plotting routines
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path

%% Grid and fracture lines

celldim = [90 90];
physdim = [9 9];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

fl = [20, 45, 70, 45; 
      45, 20, 45, 70]./10 ; % fractures lines in [x1 y1 x2 y2] format.

%% Process fracture lines

dispif(mrstVerbose, 'Processing user input...\n\n');
[G,fracture] = processFracture2D(G,fl);
fracture.aperture = 1/25; % Fracture aperture
figure;
plotFractureLines(G,fracture,'lines');
axis tight; 
box on

%% Compute CI and construct fracture grid

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 0.05; cell_size = 0.05; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); box on

%% Set rock properties in fracture and matrix

dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
G.rock.perm = ones(G.cells.num,1) * darcy;
K_frac = 10000; % Darcy
G = makeRockFrac(G, K_frac);

%% Define fluid properties

fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1000*kilogram/meter^3);

%% Define fracture connections as NNC and compute the transmissibilities

[G,T] = defineNNCandTrans(G,F,fracture);

%% Initialize state variables

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 5*barsa);

% Get A matrix without source
[A,q] = getSystemIncompTPFA(state, G, T, fluid, 'use_trans', true);

%% Setup multiscale grids

dispif(mrstVerbose, 'Defining coarse grids and interaction regions...\n\n');

coarseDims = [9 9]; 
dof_frac = 4; % Number of coarse blocks in the fracture grid
[CG, CGf] = getRsbGridsHFM(G, fracture.network, 'coarseDims', coarseDims,...
            'dof_frac',dof_frac);

clf; plotFractureCoarseGrid2D(G,CG.partition,F)

%% Add BC

bc = [];
bc  = pside(bc, G, 'LEFT', 10*barsa);
bc  = pside(bc, G, 'RIGHT', 1*barsa);

%% Incompressible 1-phase FS

W = []; % 
dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state_fs = incompTPFA(state, G, T, fluid,  ...
            'bc',bc, 'MatrixOutput', true, 'use_trans',true);

%% Compute basis functions

dispif(mrstVerbose, 'Computing basis functions...\n\n');
basis_sb = getMultiscaleBasis(CG, A, 'type', 'rsb');
clf; plotToolbar(G,basis_sb.B); 
axis tight; c = colormap(jet);  
c(1,:) = [1 1 1]; colormap(c); colorbar; 
title('Basis Functions in the matrix');

%% Compute multiscale solution

dispif(mrstVerbose, 'Computing multiscale solution...\n\n');
[state_ms,~] = incompMultiscale(state, CG, T, fluid, basis_sb,...
               'bc', bc,'use_trans',true);

%% Solve using MS-ILU and MS-GMRES

fn = getSmootherFunction('type', 'ilu');

[~,report] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'bc', bc, 'use_trans',true, 'tolerance', 1e-8, 'iterations', 100,...
    'useGMRES', false, 'reconstruct', true, 'getSmoother', fn);

[~,report2] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'bc', bc, 'use_trans',true, 'tolerance', 1e-8, 'iterations', 100,...
    'useGMRES', true, 'reconstruct', true, 'getSmoother', fn);

%% Plot results

figure;
plotToolbar(G, state_fs.pressure)
colormap jet
view(90, 90); colorbar
axis tight off
title('Fine scale')

figure;
plotToolbar(G, state_ms.pressure)
colormap jet
view(90, 90); colorbar
axis tight off
title('F-MsRSB')

L1 = abs(state_ms.pressure-state_fs.pressure)./state_fs.pressure;
figure;
plotToolbar(G, L1)
colormap jet
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