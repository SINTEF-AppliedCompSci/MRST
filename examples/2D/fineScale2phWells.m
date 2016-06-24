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

celldim = [50 20];
physdim = [500 200];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

fl = [80,  100, 270, 180;...
      130, 160, 340,  40;...
      260,  40, 420, 150] ; % fractures lines in [x1 y1 x2 y2] format.

%% Process fracture lines

dispif(mrstVerbose, 'Processing user input...\n\n');
[G,fracture] = processFracture2D(G, fl, 'verbose', mrstVerbose);
fracture.aperture = 1/25; % Fracture aperture
figure;
plotFractureLines(G,fracture,'lines');
axis equal tight; 
box on

%% Compute CI and construct fracture grid

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 5; cell_size = 10; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); 
axis equal; box on

%% Set rock properties in fracture and matrix

dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
p = gaussianField(celldim, [0.2 0.4], [11 3], 2.5);
K = p.^3.*(1e-5)^2./(0.81*72*(1-p).^2);
G.rock.poro = p(:);
G.rock.perm = K(:);
K_frac = 10000; % Darcy
G = makeRockFrac(G, K_frac, 'porosity', 0.8);
clf; plotToolbar(G, G.rock);

%% Define fluid properties

fluid = initSimpleFluid('mu' , [   1,  1] .* centi*poise     , ...
    'rho', [1000, 700] .* kilogram/meter^3, ...
    'n'  , [   2,   2]);

%% Define fracture connections as NNC and compute the transmissibilities

[G,T] = defineNNCandTrans(G,F,fracture);

%% Initialize state variables

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 0);
state.wellSol = initWellSol(G, 0);

%% Add wells

inj = 1;
prod = celldim(1)*celldim(2);
wellRadius = 0.1;

W = addWell([], G.Matrix, G.Matrix.rock, inj,'InnerProduct', 'ip_tpf','Type', ...
    'rate', 'Val', 10, 'Radius', wellRadius, 'Comp_i', [1, 0]);
W = addWell(W, G.Matrix, G.Matrix.rock, prod, 'InnerProduct', 'ip_tpf', 'Type', ...
    'bhp' , 'Val', 50*barsa, 'Radius', wellRadius, 'Comp_i', [0, 1]);

%% Compute initial pressure

dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state = incompTPFA(state, G, T, fluid, 'Wells', W, 'MatrixOutput', true, 'use_trans',true);

%% Plot results

figure;
plotToolbar(G, state.pressure)
colormap jet
view(0, 90); colorbar
title('Initial pressure');

%% Incompressible Two-Phase Flow

pv     = poreVolume(G,G.rock);
nt     = 60;
Time   = 0.5*(sum(pv)/state.wellSol(1).flux);
dT     = Time/nt;
dTplot = Time/3;
N      = fix(Time/dTplot);

pvi = zeros(nt,1);
sol = cell(nt,1);

t  = 0;
count = 1;
while t < Time,
    state = implicitTransport(state, G, dT, G.rock, fluid, 'wells', W, 'Trans', T,'verbose',true);
    
    % Check for inconsistent saturations
    assert(max(state.s(:,1)) < 1+eps && min(state.s(:,1)) > -eps);

    % Update solution of pressure equation.
    state  = incompTPFA(state, G, T, fluid, 'wells', W, 'use_trans',true);
    sol{count,1} = state;
    
    t = t + dT;
    pvi(count) = 100*(sum(state.wellSol(1).flux)*t)/sum(pv);
    count = count + 1;
    
end

%% Plot saturations

figure; plotToolbar(G,sol); 
colormap(flipud(gray)); caxis([0 1]); axis tight equal