%{
Two-phase example with a horizontal producer and injector modeling water
injection in a 3-dimensional fractured porous media using the HFM module.
The 3D solvers do not have the capability of handling intersecting fracture
planes.
%}

mrstModule add hfm;             % hybrid fracture module
mrstModule add coarsegrid;      % functionality for coarse grids
mrstModule add new-multiscale;  % MsRSB solvers
mrstModule add mrst-gui;        % plotting routines

%% Grid and fracture(s)

celldim = [25, 25, 25];
physdim = [100 100 100];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);
GlobTri = globalTriangulation(G);
fracplanes = struct;
fracplanes(1).points = [50 25 25;
                        50 25 75;
                        50 75 75;
                        50 75 25];
fracplanes(1).aperture = 1/50;
checkIfCoplanar(fracplanes)

%% Process fracture(s)

dispif(mrstVerbose, 'Processing user input...\n\n');
[G,fracplanes] = preProcessingFractures(G, fracplanes, ...
                 'GlobTri', GlobTri, ...
                 'fractureCellSize', 0.3);
fprintf('\nProcessing complete...\n');

figure; plotGrid(G,'facealpha',0);
for i = 1:numel(fieldnames(G.FracGrid))
    plotGrid(G.FracGrid.(['Frac',num2str(i)]));
end
view(15,20);

%% Set rock properties in fracture and matrix

dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');

G.rock.perm = ones(G.cells.num,1)*darcy();
G.rock.poro = 0.2*ones(G.cells.num, 1);
K_frac = 10000; % Darcy
poro_frac = 0.5;
G = makeRockFrac(G, K_frac, 'permtype','homogeneous','porosity',poro_frac);

%% Define fluid properties

fluid = initSimpleFluid('mu' , [   1,  1] .* centi*poise     , ...
    'rho', [1000, 700] .* kilogram/meter^3, ...
    'n'  , [   2,   2]);

%% Assemble global grid and compute transmissibilities

G = assembleGlobalGrid(G);
G = computeEffectiveTrans(G);
G.faces.tag = zeros(G.faces.num,1);
x = ismember(G.faces.neighbors,G.nnc.cells,'rows');
G.faces.tag(x) = 1; % Tag NNC faces. 
%
T = computeTrans(G, G.rock);
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
T = [T;G.nnc.T];

%% Initialize state variables

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 0);
state.wellSol = initWellSol(G, 0);
[A,q] = getSystemIncompTPFA(state, G, T, fluid, 'use_trans', true); 

%% Add wells

[nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
cellinj = nx*ny*(nz-1)+1:nx:nx*ny*nz;
cellprod = nx:nx:nx*ny;
W   = addWell([], G.Matrix, G.Matrix.rock, flipud(cellinj), 'Type', 'rate',...
    'Val', 500/day, 'Sign',1, 'Comp_i', [1, 0], 'Name', 'Injector');
W   = addWell(W, G.Matrix, G.Matrix.rock, cellprod, 'Type', 'bhp', ...
    'Val', 100*barsa, 'Sign', -1, 'Comp_i', [0, 1], 'Name', 'Producer');

%% Setup multiscale grids

nw = struct; 
for i = 1:numel(fracplanes)
    nw(i).lines = i;
end
% Partition matrix

coarseDims = [5 5 5];
pm = partitionMatrix(G, 'coarseDims', coarseDims, 'use_metis', false);
CGm = getRsbGridsMatrix(G, pm, 'Wells', W);

% Partition fracture

coarseDimsF = [3 4];
p  = partitionFracture(G, pm, nw, 'partition_frac'   , true   , ...
    'use_metisF'       , false  , ...
    'coarseDimsF'      , coarseDimsF );

% p = processPartition(G,compressPartition(p));
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

%% Plot coarsegrid

clf; % plotToolbar(G.Matrix,G.Matrix.rock);
colormap(jet); view(15,20);
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
clf; plotToolbar(G,basis_sb.B); view(-135,30)
axis tight; c = colormap(jet);
c(1,:) = [1 1 1]; colormap(c); colorbar;
title('Basis Functions in the matrix');

%% Compute multiscale solution

dispif(mrstVerbose, 'Computing multiscale solution...\n\n');
[state_ms,~] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'Wells', W,'use_trans',true);

%% Plot initial pressure

figure;
plotToolbar(G, state_fs.pressure)
colormap jet; colorbar
view(15,20)
axis tight off
title('Fine scale')

figure;
plotToolbar(G, state_ms.pressure)
colormap jet; colorbar
view(15,20)
axis tight off
title('F-MsRSB')

L1 = abs(state_ms.pressure-state_fs.pressure)./state_fs.pressure;
figure;
plotToolbar(G, L1)
plotWell(G,W); view(15,20);
colormap jet; colorbar
axis tight off
L1_eq = '$$ \frac{| P_i^{fs}-P_i^{f-msrsb} | }{ P_i^{fs} } $$';
title(L1_eq,'interpreter','latex');

%% Incompressible Two-Phase Flow

pv     = poreVolume(G,G.rock);
nt     = 30;
t90   = 0.9*(sum(pv)/abs(sum(state_fs.wellSol(1).flux)));
Time   = t90;
dT     = Time/nt;

pvi = zeros(nt,1);
sol_fs = cell(nt,1); sol_ms = cell(nt,1);
e = zeros(nt,1);

t  = 0;
B = basis_sb.B;
R = controlVolumeRestriction(CG.partition);
count = 1;
while t < Time,
    state_fs = implicitTransport(state_fs, G, dT, G.rock, fluid, 'wells', W, 'Trans', T,'verbose',true);
    state_ms = implicitTransport(state_ms, G, dT, G.rock, fluid, 'wells', W, 'Trans', T);
    % Check for inconsistent saturations
    s = [state_fs.s(:,1); state_ms.s(:,1)];
    assert(max(s) < 1+eps && min(s) > -eps);

    % Update solution of pressure equation.
    state_fs  = incompTPFA(state_fs , G, T, fluid, 'wells', W, 'use_trans',true);

    %-------------------------------Multiscale----------------------------%
    A = getSystemIncompTPFA(state_ms, G, T, fluid, 'use_trans', true);
    B = iteratedJacobiBasis(A, CG, 'interpolator', basis_sb.B);
    basis_sb = struct('B', B, 'R', R);
    state_ms = incompMultiscale(state_ms, CG, T, fluid, basis_sb, 'Wells', W,...
        'use_trans',true);
    %---------------------------------------------------------------------%
    
    sol_fs{count,1} = state_fs; sol_ms{count,1} = state_ms;
    % Increase time
    t = t + dT;
    pvi(count) = 100*(sum(state_fs.wellSol(1).flux)*t)/sum(pv);
    e(count,1) = sum(abs(state_fs.s(:,1) - state_ms.s(:,1)).*pv)/sum(pv.*state_fs.s(:,1));
    
    fprintf([num2str(pvi(count)), '%% PV injected \n']);
    count = count + 1;
end

%% Plot saturations

close all;
figure; plotToolbar(G,sol_fs); colormap(flipud(gray)); caxis([0 1]);view(15,20)
figure; plotToolbar(G,sol_ms); colormap(flipud(gray)); caxis([0 1]);view(15,20)

%% Plot error in saturation 

figure;
plot(pvi,e*100, '--+b');
ylabel('e [%]')
xlabel('PVI [%]'); 
set(gca,'FontSize',18,'XGrid','on','YGrid','on');
axis tight

e_eq = '$$ e = \frac{ \sum ( |S_w^{fs}-S_w^{f-msrsb}| \times pv) }{ \sum (S_w^{fs} \times pv) } $$';
title(e_eq,'interpreter','latex');