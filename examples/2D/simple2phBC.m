%{
Two-phase 2D example with water injection from the left boundary of a
rectangular domain containing 2 intersecting fractures.
%}

close all;

%% Grid and fracture lines

celldim = [99 99];
physdim = [99 99];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

fl = [ 22    22    77    77; 22    77    77    22];

%% Process fracture lines

dispif(mrstVerbose, 'Processing user input...\n\n');
[G,fracture] = processFracture2D(G,fl);
fracture.aperture = 1/25; % Fracture aperture
figure;
plotFractureLines(G,fracture,'lines');
box on

%% Compute CI and construct fracture grid

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 0.1; cell_size = 0.5; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); box on

%% Set rock properties in fracture and matrix

dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
G.rock.perm = ones(G.cells.num,1)*darcy;
G.rock.poro = 0.2*ones(G.cells.num, 1);
K_frac = 10000; % Darcy
poro_frac = 0.5;
G = makeRockFrac(G, K_frac, 'permtype','homogeneous','porosity',poro_frac);

%% Define fluid properties

fluid = initSimpleFluid('mu' , [   1,  2] .* centi*poise     , ...
    'rho', [1000, 700] .* kilogram/meter^3, ...
    'n'  , [   2,   2]);

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

dispif(mrstVerbose, 'Computing basis functions...\n\n');
basis_sb = getMultiscaleBasis(CG, A, 'type', 'rsb');
clf; plotToolbar(G,basis_sb.B);
axis tight; c = colormap(jet);
c(1,:) = [1 1 1]; colormap(c); colorbar;
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

pv     = poreVolume(G,G.rock);
nt     = 30;
t90    = 2*(sum(pv)/sum(state_fs.flux(left)));
Time   = t90;
dT     = Time/nt;
dTplot = Time/3;
N      = fix(Time/dTplot);

state_ms.rhs = q;
t  = 0; plotNo = 1; hfs = 'Reference Saturation: '; hms = 'F-MsRSB Saturation: ';
dispif(mrstVerbose, 'Solving Transport...\n\n');
close all
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
