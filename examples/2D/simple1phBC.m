close all; 

%% Grid and fracture lines

celldim = [225 225];
physdim = [9 9];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

fl = [2,4.5,7,4.5;4.5,2,4.5,7];

%% Process fracture lines into grid

dispif(mrstVerbose, 'Processing user input...\n\n');
[G,fracture] = processFracture2D(G,fl);
fracture.aperture = 1/25; % Fracture aperture
figure; 
plotFractureLines(G,fracture,'lines');
box on

%% Compute CI and construct fracture grid

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 0.1; cell_size = 0.2; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); box on

%% Set rock properties in fracture and matrix

dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
G.rock.perm = ones(G.cells.num,1)*darcy();
K_frac = 10000; % Darcy
G = makeRockFrac(G, K_frac);
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1000*kilogram/meter^3);

%% Define fracture connections as NNC and compute the transmissibilities

[G,T] = defineNNCandTrans(G,F,fracture);

%% Initialize state variables

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 5*barsa());

% Get A matrix without source
A = incompTPFA(state, G, T, fluid, 'MatrixOutput', true, 'use_trans',true); 
A = A.A; A(1,1) = A(1,1)/2; % undo magic done in incompTPFA


%% Setup Multiscale Grids

dispif(mrstVerbose, 'Defining coarse grids and interaction regions...\n\n');

coarseDims = [9 9]; 
dof_frac = 3; % Number of coarse blocks in the fracture grid
[CG, CGf] = getRsbGridsHFM(G, fracture.network, 'coarseDims', coarseDims,...
            'dof_frac',dof_frac);

clf; plotFractureCoarseGrid2D(G,CG.partition,F)

%% Add BC

bc = [];
bc  = pside(bc, G, 'LEFT', 10*barsa());
bc  = pside(bc, G, 'RIGHT', 1*barsa());

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

%% Compute Multiscale solution

dispif(mrstVerbose, 'Computing multiscale solution...\n\n');
[state_ms,~] = incompMultiscale(state, CG, T, fluid, basis_sb,...
               'bc', bc,'use_trans',true);

%% Solve using MS-GMRES
fn = getSmootherFunction('type', 'ilu');

[~,report] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'bc', bc, 'use_trans',true, 'tolerance', 1e-12, 'iterations', 100,...
    'useGMRES', true, 'reconstruct', true, 'getSmoother', fn);

%% Plot results

figure;
plotToolbar(G, state_fs.pressure)
colormap jet
view(90, 90)
axis tight off
title('Fine scale')

figure;
plotToolbar(G, state_ms.pressure)
colormap jet
view(90, 90)
axis tight off
title('F-MsRSB')

%% Plot Convergence
figure;
semilogy(report.resvec, '*-')
title('Convergence of GMRES')