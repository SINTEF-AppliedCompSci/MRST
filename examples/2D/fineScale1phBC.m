%{
Single-phase 2D example with Dirichlet boundary conditions and a horizontal
central fracture in the centre comparing the embedded discrete fracture
model to a fully resolved simulation where the fracture and matrix grid
blocks are of the same size.
%}

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path

%% Grid and fracture lines

celldim = [90 90];
physdim = [9 9];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

fl = [20, 45, 70, 45]./10 ; % fractures lines in [x1 y1 x2 y2] format.

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
min_size = 0.05; cell_size = 0.1; % minimum and average cell size.
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

%% Define completely resolved grid

celldim = [225 225];
physdim = [9 9];
Gr = cartGrid(celldim, physdim);
Gr = computeGeometry(Gr);
Gr.rock.perm = ones(Gr.cells.num,1) * darcy;
fractureCells = 225*112+50:225*112+175;
Gr.rock.perm(fractureCells) = 10000*darcy;
clf; plotCellData(Gr,Gr.rock.perm);

%% Initialize state variables

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 5*barsa);
state_r  = initResSol (Gr, 5*barsa);

%% Add BC

bc = [];
bc  = pside(bc, G, 'LEFT', 10*barsa);
bc  = pside(bc, G, 'RIGHT', 1*barsa);

bcr = [];
bcr  = pside(bcr, Gr, 'LEFT', 10*barsa);
bcr  = pside(bcr, Gr, 'RIGHT', 1*barsa);

%% Compute pressure

dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state = incompTPFA(state, G, T, fluid, 'bc', bc, 'use_trans',true);
state_r = incompTPFA(state_r, Gr, computeTrans(Gr,Gr.rock), fluid, 'bc', bcr);

%% Plot results

figure;
plotToolbar(G, state.pressure)
colormap jet
view(0, 90); colorbar
title('Embedded fracture model');

figure;
plotToolbar(Gr, state_r.pressure)
colormap jet
view(0, 90); colorbar
title('Fully resolved solution')