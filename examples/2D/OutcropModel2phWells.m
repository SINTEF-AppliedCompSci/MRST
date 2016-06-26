%{
2ph example modeling water injection into a highly fractured oil-filled
reservoir. The fracture-network has been extracted from an outcrop model.

K. Bisdom, B. D. M. Gauthier, G. Bertotti, N. J. Hardebol. Calibrating
discrete fracture-network models with a carbonate three-dimensional outcrop
fracture network: Implications for naturally fractured reservoir modeling.
AAPG Bulletin, 24 (2014) 1351-1376.

Notice that you need to have Metis installed to get this example to work.
To get Metis working, you also need to set the global variable METISPATH.
This can be done in your 'startup_user.m' file.
%}

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add coarsegrid;      % functionality for coarse grids
mrstModule add ad-core;         % NNC support for coarse grids
mrstModule add new-multiscale;  % MsRSB solvers
mrstModule add mrst-gui;        % plotting routines
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path

%% Grid and fracture lines

celldim = [100 100];
physdim = [1000 1000];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

load examples/data/brazil_fractures

%% Process fracture lines

% dispif(mrstVerbose, 'Processing user input...\n\n'); tic;
% [G,fracture] = processFracture2D(G,fl); toc
% fracture.aperture = 1/25; % Fracture aperture

% The processing algorithm might take time. Comment the above 3 lines and
% load the files brazil_fractures_processed and brazil_grid_processed

load examples/data/brazil_fractures_processed
load examples/data/brazil_grid_processed

figure;
plotFractureLines(G,fracture);
box on

%% Compute CI and construct fracture grid

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 5; cell_size = 10; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); box on

%% Set rock properties in fracture and matrix

dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
load examples/data/brazil_perm

G.rock = makeRock(G,p(:),0.2);
K_frac = 1000; % Darcy
poro_frac = 0.5;
for i = 1:numel(fieldnames(G.FracGrid))
    G.FracGrid.(['Frac',num2str(i)]).rock.perm = darcy*K_frac*ones(G.FracGrid.(['Frac',num2str(i)]).cells.num,1);
    G.FracGrid.(['Frac',num2str(i)]).rock.poro = poro_frac*ones(G.FracGrid.(['Frac',num2str(i)]).cells.num,1);
end
clf; plotToolbar(G,G.rock); colormap(jet); colorbar

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
[A,q] = getSystemIncompTPFA(state, G, T, fluid, 'use_trans', true); 

%% Setup multiscale grids - This step may take some time

dispif(mrstVerbose, 'Defining coarse grids and interaction regions...\n\n');
coarseDims = [15 15]; % coarsening factor in each direction

nw = fracture.network;
numfc = zeros(1,numel(fieldnames(G.FracGrid)));
numnc = zeros(1,numel(nw));
for i = 1:numel(nw)
    for j = 1:numel(nw(i).lines)
        numfc(nw(i).lines(j)) = G.FracGrid.(['Frac',num2str(nw(i).lines(j))]).cells.num;
    end
    numnc(i) = sum(numfc(nw(i).lines));
end

dof_frac = ceil(numnc./20); tic; % Fracture dof per fracture network
[CG, CGf] = getRsbGridsHFM(G, nw, 'coarseDims', coarseDims,...
    'dof_frac',dof_frac); 

dispif(mrstVerbose, '\nMultiscale grids defined...\n\n');
toc
clf; plotFractureCoarseGrid2D(G,CG.partition,F)

%% Add wells

inj = celldim(2)*(celldim(1)-1)+1;
prod = celldim(1);
wellRadius = 0.1;

W = addWell([], G.Matrix, G.Matrix.rock, inj,'InnerProduct', 'ip_tpf','Type', ...
    'rate', 'Val', 100/day,'Radius', wellRadius, 'Comp_i', [1, 0]);
W = addWell(W, G.Matrix, G.Matrix.rock, prod, 'InnerProduct', 'ip_tpf', 'Type', ...
    'bhp' , 'Val', 10*barsa, 'Radius', wellRadius, 'Comp_i', [0, 1]);

%% Compute initial pressure

dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state_fs = incompTPFA(state, G, T, fluid,  ...
    'Wells',W, 'MatrixOutput', true, 'use_trans',true);

%% Compute basis functions

dispif(mrstVerbose, 'Computing basis functions...\n\n');
basis_sb = getMultiscaleBasis(CG, A, 'type', 'rsb');
clf; plotToolbar(G,basis_sb.B);
axis tight; c = colormap(jet);
colormap(c); colorbar;
title('Basis functions plotted in the matrix');

%% Compute multiscale solution

dispif(mrstVerbose, 'Computing multiscale solution...\n\n');
[state_ms,~] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'Wells', W,'use_trans',true);

%% Incompressible Two-Phase Flow

pv     = poreVolume(G,G.rock);
nt     = 90;
Time   = 0.75*(sum(pv)/state_fs.wellSol(1).flux);
dT     = Time/nt;
dTplot = Time/3; 
N      = fix(Time/dTplot);

pvi = zeros(nt,1);
sol_fs = cell(nt,1); sol_ms = cell(nt,1);
e = zeros(nt,1); pms = zeros(nt,3); pfs = zeros(nt,3);

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
    B = iteratedJacobiBasis(A, CG, 'interpolator', B); 
    basis_sb = struct('B', B, 'R', R);
    state_ms = incompMultiscale(state_ms, CG, T, fluid, basis_sb, 'Wells', W,'use_trans',true);

    %---------------------------------------------------------------------%
    
    sol_fs{count,1} = state_fs; sol_ms{count,1} = state_ms;
    % Increase time
    t = t + dT;
    pvi(count) = 100*(sum(state_fs.wellSol(1).flux)*t)/sum(pv);
    
    e(count,1) = sum(abs(state_fs.s(:,1) - state_ms.s(:,1)).*pv)/sum(pv.*state_fs.s(:,1));
    pfs(count,1) = state_fs.s(W(2).cells,1);
    pms(count,1) = state_ms.s(W(2).cells,1);
    
    count = count + 1;
    
end

%% Plot saturations

plotNo = 1;
figure; hold on
boundary = any(G.Matrix.faces.neighbors==0,2);
facelist = 1:G.Matrix.faces.num;
bfaces = facelist(boundary);
for i = nt/3:nt/3:nt
    state_fs = sol_fs{i,1}; state_ms = sol_ms{i,1};
    heading = [num2str(round(pvi(i),1)),  ' % PVI'];
    % Plot saturation
    r = 0.01;
    subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.48]), cla
    plotToolbar(G,state_fs.s(:,1));
    for j = 1:size(fl,1)
        line(fl(j,1:2:3),fl(j,2:2:4),'Color','r','LineWidth',0.5);
    end
    plotFaces(G,bfaces,'k','linewidth',1)
    axis square off, 
    title(['Reference Saturation: ', heading]);
    colormap(flipud(gray)), view(0,90); caxis([0 1]);
    
    subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.48]), cla
    plotToolbar(G,state_ms.s(:,1));
    for j = 1:size(fl,1)
        line(fl(j,1:2:3),fl(j,2:2:4),'Color','r','LineWidth',0.5);
    end
    plotFaces(G,bfaces,'k','linewidth',1)
    axis square off, 
    title(['F-MsRSB Saturation: ',  heading]);
    colormap(flipud(gray)), view(0,90); caxis([0 1]);
    plotNo = plotNo+1;
    
end

%% Plot water saturation at producer 

figure;
plot(pvi,pfs(:,1),'-o',pvi,pms(:,1),'--*');
leg = legend('Fine-scale','Multiscale','Location','Best');
ylabel('Saturation at producer');
xlabel('PVI [%]'); 
set(gca,'XGrid','on','YGrid','on');
axis tight

%% Plot error in saturation 

figure;
plot(pvi,e*100, '--+b');
ylabel('e [%]')
xlabel('PVI [%]'); 
set(gca,'XGrid','on','YGrid','on');
axis tight

e_eq = '$$ e = \frac{ \sum ( |S_w^{fs}-S_w^{f-msrsb}| \times pv) }{ \sum (S_w^{fs} \times pv) } $$';
title(e_eq,'interpreter','latex');