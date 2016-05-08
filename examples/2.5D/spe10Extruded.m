%{
2ph example with a statistical fracture distribution. Matrix permeability
is sampled from the SPE10 test case onto a lower resolution grid.

K. Bisdom, B. D. M. Gauthier, G. Bertotti, N. J. Hardebol. Calibrating
discrete fracture-network models with a carbonate three-dimensional outcrop
fracture network: Implications for naturally fractured reservoir modeling.
AAPG Bulletin, 24 (2014) 1351-1376.

%}

close all;
mrstModule add spe10

%% Grid and fracture lines

celldim = [30 110];
G = cartGrid(celldim);
G = computeGeometry(G);
load datasets/statistical_fractures

layers = 40;
flayers = 11:30;

%% Process fracture lines

dispif(mrstVerbose, 'Processing user input...\n\n');
a = 1/25;
[G,fracture] = processFracture2D(G,fl); fracture.aperture = a;
figure;
plotFractureLines(G,fracture,'lines');
axis tight; box on

%% Compute CI and construct fracture grid

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 0.5; cell_size = 1; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); box on

%% Make layered grid

Gl = makeLayers(G,layers,flayers);

%% Set rock properties in fracture and matrix

dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
rock = SPE10_rock(1:85);
perm = sampleFromBox(Gl,reshape(rock.perm(:,1).*milli.*darcy,60,220,[]));
poro = sampleFromBox(Gl,reshape(rock.poro,60,220,[]));
% Since we don't use actnum to denote N/G or active cells we ignore cells
% with very low porosity.
poro(poro<0.01) = 0.01;

Gl.rock = makeRock(Gl, perm(:), poro(:));
K_frac = 1000; % Darcy
poro_frac = 0.8;
Gl = makeRockFrac(Gl, K_frac, 'permtype','homogeneous','porosity',poro_frac);

clf; plotToolbar(Gl,Gl.rock); 
colormap(jet); view (-135,30);
axis tight equal off

%% Define fluid properties

fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);

%% Define fracture connections as NNC and compute the transmissibilities

[Gl,T] = makeNNCextruded(G,Gl,F,fracture,flayers);
G = Gl; clear Gl

%% Initialize state variables

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 0);
state.wellSol = initWellSol(G, 0);
[A,q] = getSystemIncompTPFA(state, G, T, fluid, 'use_trans', true);

%% Set up 5 spot

W = [];
% Add an injector in the middle of the domain based on the midpoints of the
% cartesian indices
midpoint = round(G.cartDims ./ 2);

% Central injector
W = verticalWell(W, G.Matrix, G.Matrix.rock, midpoint(1), midpoint(2), [], ...
    'Type', 'rate',...  
    'Val', 1000*stb/day,...        
    'Radius', .125*meter, ...    
    'InnerProduct', 'ip_tpf', ...
    'Name', 'I1', ...
    'Comp_i', [1 0]);           

% Producer at each corner
h_ind = [1, G.cartDims(1)];
v_ind = [1, G.cartDims(2)];
count = 1;
for i = 1:2
    for j = 1:2
        W = verticalWell(W, G.Matrix, G.Matrix.rock, h_ind(i), v_ind(j), [], ...
            'Type', 'bhp',...
            'Val', 1000*psia, ...
            'Radius', .125*meter, ...
            'InnerProduct', 'ip_tpf',...
            'Name', ['P', num2str(count)], ...
            'Comp_i', [0, 1]);
        count = count + 1;
    end
end

%% Setup multiscale grids

dispif(mrstVerbose, 'Defining coarse grids and interaction regions...\n\n');
G.type{1,1} = 'layered';
nw = fracture.network;
numfc = zeros(1,numel(fieldnames(G.FracGrid)));
numnc = zeros(1,numel(nw));
for i = 1:numel(nw)
    for j = 1:numel(nw(i).lines)
        numfc(nw(i).lines(j)) = G.FracGrid.(['Frac',num2str(nw(i).lines(j))]).cells.num;
    end
    numnc(i) = sum(numfc(nw(i).lines));
end
cf = 80; % Coarsening factor per network

dof_frac = ceil(numnc./cf);
coarseDims = [6 22 8]; tic; % coarsening factor in each direction
[CG,CGf] = getRsbGridsHFM(G, fracture.network, 'coarseDims', coarseDims, ...
           'dof_frac', dof_frac, 'sysMatrix', A, 'paddedPartition', true, ...
           'Wells', W);      
dispif(mrstVerbose, '\nMultiscale grids defined...\n\n');toc

%% Plot coarsegrid

clf; plotToolbar(G.Matrix,G.Matrix.rock);
colormap(jet); view(-135,30);
axis tight equal off
plotGrid(CG, 'facealpha', 0, 'linewidth', 2);
plotWell(G,W);

%% Compute initial pressure

dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state_fs = incompTPFA(state, G, T, fluid,  ...
    'Wells', W, 'MatrixOutput', true, 'use_trans',true);

%% Compute basis functions

dispif(mrstVerbose, 'Computing basis functions...\n\n');
basis_sb = getMultiscaleBasis(CG, A, 'type', 'rsb');
clf; plotToolbar(G,basis_sb.B);
axis tight; c = colormap(jet);
c(1,:) = [1 1 1]; colormap(c); colorbar;
title('Basis functions plotted in the matrix');

%% Compute Multiscale solution

tol = 1e-3; % choose tolerance 
iter = 2;  % choose number of iterations
fn = getSmootherFunction('type', 'ilu');

%-------------------------------F-MsRSB-----------------------------------%
[state_ms,report0] = incompMultiscale(state, CG, T, fluid, basis_sb, ...
                     'Wells', W, 'use_trans',true);
%----------------------------F-MsRSB + ILU(0)-----------------------------%
[state_i0,reporti0] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'Wells', W, 'use_trans',true, 'reconstruct', true,'getSmoother', fn, ...
    'iterations', iter, 'tolerance', tol);

%% Incompressible Two-Phase Flow

pv     = poreVolume(G,G.rock);
nt     = 10;
Time   = 10*100*day;
dT     = Time/nt;

pvi = zeros(nt,1);
sol_fs = cell(nt,1); sol_ms = cell(nt,1); 
sol_i0 = cell(nt,1);
t  = 0;
B = basis_sb.B; 
basis_sbi0 = basis_sb; 
R = controlVolumeRestriction(CG.partition);
count = 1;
while t < Time,
    state_fs = implicitTransport(state_fs, G, dT, G.rock, fluid, 'wells', W, 'Trans', T, 'verbose', true);
    state_ms = implicitTransport(state_ms, G, dT, G.rock, fluid, 'wells', W, 'Trans', T);
    state_i0 = implicitTransport(state_i0, G, dT, G.rock, fluid, 'wells', W, 'Trans', T);
    % Check for inconsistent saturations
    
    s = [state_fs.s(:,1); state_ms.s(:,1)];
    assert(max(s) < 1+eps && min(s) > -eps);
    
    % Update solution of pressure equation.
    state_fs  = incompTPFA(state_fs , G, T, fluid, 'wells', W, 'use_trans',true, 'verbose', true);
    
    %-------------------------------Multiscale----------------------------%
    A = getSystemIncompTPFA(state_ms, G, T, fluid, 'use_trans', true);
    Ai0 = getSystemIncompTPFA(state_i0, G, T, fluid, 'use_trans', true);
    
    B = iteratedJacobiBasis(A, CG, 'interpolator', basis_sb.B);
    basis_sb = struct('B', B, 'R', R);
    
    Bi0 = iteratedJacobiBasis(Ai0, CG, 'interpolator', basis_sbi0.B);
    basis_sbi0 = struct('B', Bi0, 'R', R);
    
    state_ms = incompMultiscale(state_ms, CG, T, fluid, basis_sb, 'Wells', W,...
        'use_trans',true);
   
    state_i0 = incompMultiscale(state_i0, CG, T, fluid, basis_sbi0,...
        'Wells', W, 'use_trans',true, 'reconstruct', true, 'getSmoother', fn, ...
        'iterations', iter, 'tolerance', tol);
    %---------------------------------------------------------------------%
    
    sol_fs{count,1} = state_fs; sol_ms{count,1} = state_ms;
    sol_i0{count,1} = state_i0; 
    
    % Increase time
    t = t + dT;
    pvi(count) = 100*(sum(state_fs.wellSol(1).flux)*t)/sum(pv);
    fprintf([num2str(pvi(count)), '%% PV injected \n']);
    count = count + 1;
end

%% Plot saturations

close all;
figure; plotToolbar(G,sol_fs); colormap(flipud(gray)); caxis([0 1]);
figure; plotToolbar(G,sol_ms); colormap(flipud(gray)); caxis([0 1]);

%% Compute well solutions

wellSol_fs = convertIncompWellSols(W, sol_fs);
wellSol_ms = convertIncompWellSols(W, sol_ms);
wellSol_i0 = convertIncompWellSols(W, sol_i0);
for i = 1:numel(wellSol_fs)
wellSol_fs{i}(1).watercut = 1;
wellSol_fs{i}(2).watercut = wellSol_fs{i}(2).qWs/wellSol_fs{i}(2).qOs;
wellSol_ms{i}(1).watercut = 1;
wellSol_ms{i}(2).watercut = wellSol_ms{i}(2).qWs/wellSol_ms{i}(2).qOs;
wellSol_i0{i}(1).watercut = 1;
wellSol_i0{i}(2).watercut = wellSol_i0{i}(2).qWs/wellSol_i0{i}(2).qOs;
end

%% Plot well solutions

tol = reporti0.finalResidual;
names = {'Reference','F-MsRSB', ['Tolerance = ',num2str(tol),'; Iteration(s) = ',num2str(iter)]};
linestyles = {'-',':+',':',':',':'};
plotWellSols({wellSol_fs, wellSol_ms,wellSol_i0},pvi*day,...
    'linestyles',linestyles,'datasetnames',names);
ax = gca; ax.FontSize = 24; 
ax.XColor = 'k';
ax.YColor = 'k';
ax.XLabel.String = 'PVI [%]';
