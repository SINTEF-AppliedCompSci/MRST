% simpleIncomp2D_1ph
close all
mrstModule add hfm coarsegrid incomp new-multiscale
hasLSI = checkLineSegmentIntersect();

%% Define a grid and fracture lines
% Fracture lines are to be represented by their endpoints in the form
% [x1 y1 xy y2]. Hence, for n_f fractures, the user must provide a n_f-by-4
% vector with fracture end points.

G = cartGrid([50 50]); % A 50-by-50 grid in a 20-by-20 domain
G = computeGeometry(G);
fl = [10 10 40 40; 10 40 40 10];

%% Define a fracture aperture and fracture grid size
% The user can specify a single fracture aperture for every fracture line
% or one value per fracture line. Additionally for gridding the 2D
% fracture, the user can specify a minimum and average desired cell size
% for a fracture grid cell. If a fracture line is smaller than the minimum
% grid size then the solver will override the minimum size constraints and
% assign 1 fine cell to the respective line. For more options see
% assembleFracNodes2D.

a = 1/25; % fracture aperture
min_size = 1; % minimum fracture cell size
cell_size = 2; % desired mean fracutre cell size
[G,fracture] = processFracture2D(G,fl,'verbose',true); fracture.aperture = a;

G = CIcalculator2D(G,fracture);
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,...
    'cell_size',cell_size,'verbose',true);

% see fracture grid
plotFractureNodes2D(G, F, fracture);

%% Define rock and fluid properties
% For setting fracture permeability significantly higher than matrix
% permeability, a ratio K_star = K_frac/K_matrix can be supplied. For
% heterogeneous media, this relation will not hold as the solver will
% simply generate a heterogeneous permeability set for the fractures and
% multiply it by K_star.

fprintf('Initializing rock and fluid properties...\n\n');
G.rock.perm = ones(G.cells.num,1)*darcy();
G.rock.poro = 0.2*ones(G.cells.num, 1);
K_star = 10000; % Scaling factor = K_frac/K_mat (Darcy)
poro_frac = 0.5;
G = makeRockFrac(G, K_star, 'permtype','homogeneous','rockporo',poro_frac);
fluid = initSingleFluid('mu' , 1*centi*poise, ...
    'rho', 1000*kilogram/meter^3);


%% Define fracture matrix connections as NNC's and compute their transmissibilities
% defineNNC_computeTrans calls frac_matrix_nnc, assembleGlobalGrid and
% frac_frac_nnc internally. Returned transmissibility matrix T contains 1
% transmissibility per face. Matrix grid is now stored in G.Matrix.

[G,T] = defineNNC_computeTrans(G,F,fracture);


%% Initialize state variables

pinit = 5*barsa();
state  = initResSol (G, pinit);
% Get A matrix without source
A = incompTPFA(state, G, T, fluid, 'MatrixOutput', true, 'use_trans',true); 
A = A.A; A(1,1) = A(1,1)/2; % undo magic done in incompTPFA due to no source


%% Setup Multiscale Grids
% Several partitioning options can be specified. see getRsbGrids_HFM for
% more details.

coarsen = [15 15]; % coarsening factor in each direction
dof_frac = 5; % fracture degrees of freedom at coarse scale
[CG, CGf] = getRsbGrids_HFM(G, F, fracture.network, 'coarsen', coarsen, ...
    'dof_frac',dof_frac,'sysMatrix',A);

%see multiscale grids
clf; plotFractureCoarseGrid(G,CG.partition,F);

%% Add wells

W   = verticalWell([], G.Matrix, G.Matrix.rock, 1, 1, [], 'Type', 'rate',...
    'Val', 1/day(), 'InnerProduct', 'ip_tpf', 'Sign',1);
W   = verticalWell(W, G.Matrix, G.Matrix.rock, G.Matrix.cartDims(1), G.Matrix.cartDims(2), [], ...
    'Type', 'bhp', 'Val', 1*barsa(), 'InnerProduct', 'ip_tpf', 'Sign', -1);

%% Solve fine-scale system

state_fs = incompTPFA(state, G, T, fluid,  ...
    'Wells', W, 'MatrixOutput', true, 'use_trans',true);

%% Solve multi-scale system

basis_sb = getMultiscaleBasis(CG, A, 'type', 'rsb');
state_ms = incompMultiscale(state, CG, T, fluid, basis_sb, 'Wells', W, ...
    'use_trans',true);

%% Plots

sol = struct('FS',state_fs,'MS',state_ms); solvers = {'FS', 'MS'};

clf; set(gcf,'Position',get(0,'screensize'));
%-------------------------------Patch-------------------------------------%
for figs = 1:2
    subplot(1,2,figs);
    state = getfield(sol,solvers{1,figs}); %#ok
    for i = 1:size(fl,1)
        line(fl(i,1:2:3),fl(i,2:2:4),'Color','k','LineWidth',1);
    end
    hold on
    plotToolbar(G.Matrix, state.pressure(1:G.Matrix.cells.num),'FaceAlpha',0.9,'EdgeAlpha',0.02)
    title(['Reservoir Pressure - ',solvers{1,figs}],'FontSize',15,'FontWeight','bold');
    colormap(jet)
    colorbar
    axis tight
end

%-----------------------------Time of Flight------------------------------%
figure('Position',get(0,'screensize'));
for figs = 1:2
    subplot(1,2,figs);
    state = getfield(sol,solvers{1,figs}); %#ok
    hold on
    D = computeTOFandTracer(state,G,G.rock,'Wells',W);
    tofs = log(abs(sum(D.tof,2)));
    m1 = max(min(tofs),mean(tofs)-2*std(tofs));
    m2 = min(max(tofs),mean(tofs)+std(tofs));
    plotToolbar(G.Matrix,tofs(1:G.Matrix.cells.num),'FaceAlpha',0.9,'EdgeAlpha',0.02);
    colorbar;
    colormap(jet);
    caxis([m1 m2]);
    axis tight
    for i = 1:size(fl,1)
        line(fl(i,1:2:3),fl(i,2:2:4),'Color','k','LineWidth',1);
    end
    title(['Total TOF - ',solvers{1,figs}],'FontSize',15,'FontWeight','bold'); colorbar
end