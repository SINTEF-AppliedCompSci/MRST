close all;
clearvars -except METISPATH mrstVerbose screenSize
min_size = 2; cell_size = 5;
a = 1/25; % fracture aperture
dof_frac = 1; % Coarse dof
dof_matrix = 50; % Coarse dof
fracIRtype = 4; % connectivity based
%% User defined fracture lines and appropriate grid

% nx = 100; ny = 100;
% G = cartGrid([nx ny], [100 100]);
% G = computeGeometry(G);'Frac'
% fl = [70 105 120 110; 70 85 130 105; 80 75 100 120; 115 120 115 75; 80 95 120 95; 90 120 110 105];

% load('Lines13_Network1_100by100Domain.mat');
% load('Lines10_Network5_20by20Domain.mat');
load('Hajibeygi2011_case1_GridWithFrac.mat')
% load('Hajibeygi2011_case2_GridWithFrac.mat');
% load('BrazilOutcrop_GridWithFrac.mat');
% load('Statistical5_GridWithFrac.mat');

%% Process fracture lines into grid

[G,fracture] = processFracture2D(G,fl); fracture.aperture = a;
% figure; plotFractureLines(G,fracture,'network'); % 3rd arg can be line #'s, network or none
% figure; plotMarkedCells(G, fracture,'lines'); % 3rd arg can be line #'s, network or none


%% Compute CI and grid fracture

dispif(mrstVerbose, 'Computing CI and gridding fracture...\n\n');
G = CIcalculator2D(G,fracture);
% plotCI(G,fracture);

[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
% figure; plotFractureNodes2D(G,F,fracture,show_node_num); clear show_node_num


%% Set rock permeabilities in fracture and matrix

dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
G.rock.perm = ones(G.cells.num,1)*darcy();
% G.rock.perm = randPermGen(G.cells.num);
% G.rock.perm = randi(5e10,G.cells.num,1)*darcy()./randi(5e12,G.cells.num,1);
G.rock.poro = 0.2*ones(G.cells.num, 1);
K_frac = 10000; % Scaling factor = K_frac/K_mat D
poro_frac = 0.5;
G = makeRockFrac(G, K_frac, 'permtype','homogeneous','rockporo',poro_frac);
fluid = initSingleFluid('mu' , 1*centi*poise, ...
    'rho', 1000*kilogram/meter^3);
% plotMatrixPermeability(G);


%% Define NNC and corresponding Trans

[G,T] = defineNNCandTrans(G,F,fracture);


%% Initialize state variables

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
pinit = 5*barsa(); wpinit = 10*barsa();
state  = initResSol (G, pinit);
state.wellSol = initWellSol(G, wpinit);
% Get A matrix without source
A = incompTPFA(state, G, T, fluid, 'MatrixOutput', true, 'use_trans',true); 
A = A.A; A(1,1) = A(1,1)/2; % undo magic done in incompTPFA


%% Setup Multiscale Grids

dispif(mrstVerbose, 'Defining coarse grids and interaction regions...\n\n');

coarsen = [8 8]; % coarsening factor in each direction
% set partition_frac to false if you don't want a separate interaction
% region for fractures. i.e. no rsb in fractures
[CG, CGf] = getRsbGrids_HFM(G, fracture.network, 'coarsen', coarsen, ...
    'use_metis', false,'dof_matrix',dof_matrix,'dof_frac',dof_frac,...
    'sysMatrix',A,'fracIRtype',fracIRtype,'fullyCoupled',true);

dof_matrix = max(CG.partition(1:G.Matrix.cells.num));
figure; plotFractureCoarseGrid(G,CG.partition,F)


%% Plot interaction region
% 
% IR = zeros(G.cells.num,CG.cells.num);
% for i = 1:CG.cells.num
%     IR(CG.cells.interaction{i,1},i) = 1;
% end
% figure; plotToolbar(G,IR,'EdgeAlpha',0.05); hold on
% outlineCoarseGrid(G.Matrix,CG.partition(1:G.Matrix.cells.num)); axis tight; c = colormap(flipud(winter));  
% caxis([0 1]); c(1,:) = [1 1 1]; colormap(c); 
% title('Interaction Regions','FontSize',15,'FontWeight','bold');

% IRf = zeros(CGf.parent.cells.num,CGf.cells.num);
% for i = 1:CGf.cells.num
%     IRf(CGf.cells.interaction{i,1},i) = 1;
% end
% figure; plotFracData(G, IRf, 'figTitle','Internal fracture interaction regions',...
%     'wide',true,'F',F,'CG',CGf,'plotMatrix',true,'cmap',flipud(winter));


%% Add Wells and BC

% No bc = no-flow boundaries
bc = [];
bc  = pside(bc, G, 'LEFT', 10*barsa());
bc  = pside(bc, G, 'RIGHT', 1*barsa());
radius = 1e-2;
W   = verticalWell([], G.Matrix, G.Matrix.rock, 1, 1, [], 'Type', 'rate',...
    'Val', 0.1/day(), 'InnerProduct', 'ip_tpf', 'Radius', radius, 'Sign',1);
W   = verticalWell(W, G.Matrix, G.Matrix.rock, G.Matrix.cartDims(1), G.Matrix.cartDims(2), [], ...
    'Type', 'bhp', 'Val', 0.5*barsa(), 'InnerProduct', 'ip_tpf', 'Radius', radius, 'Sign', -1);


%% Incompressible 1-phase FS

W = []; %
dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state_fs = incompTPFA(state, G, T, fluid,  ...
    'bc',bc, 'Wells', W, 'MatrixOutput', true, 'use_trans',true);
[neighborship, nnc] = getNeighbourship(G,'Topological',1);
[cellNo, cellFaces, isNNC] = getCellNoFaces(G);
nncfaces = cellFaces(isNNC);
flux_nnc = sum(state_fs.flux(nncfaces));


%% Multiscale basis and solution

dispif(mrstVerbose, 'Solving multiscale system...\n\n');

%----------------------------Basis functions------------------------------%

basis_sb = getMultiscaleBasis(CG, A, 'type', 'rsb');
figure; plotToolbar(G,basis_sb.B); axis tight; c = colormap(jet);  
c(1,:) = [1 1 1]; colormap(c); colorbar; 
title('Basis Functions','FontSize',15,'FontWeight','bold');

%-----------Dirty way to see basis functions inside fractures-------------%

figure; plotFracData(G, basis_sb.B(G.Matrix.cells.num+1:G.cells.num,1:CG.cells.num), ...
    'wide',true,'F',F,'CG',CGf,'plotMatrix',true,'cmap',jet);

%----------------------------1 MSFV operation-----------------------------%

[state_ms,report1] = incompMultiscale(state, CG, T, fluid, basis_sb, 'Wells', W, ...
    'bc', bc,'use_trans',true);

%----------------------------Use ILU Smoother-----------------------------%
fn = getSmootherFunction('type', 'ilu');

[state_ms_iter,report] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'Wells', W, 'bc', bc, 'use_trans',true, 'tolerance', 1e-6, 'iterations', 1e3,...
    'useGMRES', false, 'reconstruct', true, 'getSmoother', fn);
sol = struct('FS',state_fs,'MS',state_ms); solvers = {'FS', 'MS'};


%% Plots - pressure

dispif(mrstVerbose, 'Plots...\n\n');
% figure('Position',[screenSize(1:end-1),screenSize(end)/1.5])
% %-------------------------------Patch-------------------------------------%
% for figs = 1:2
%     subplot(1,2,figs);
%     state = getfield(sol,solvers{1,figs}); %#ok
%     for i = 1:size(fl,1)
%         line(fl(i,1:2:3),fl(i,2:2:4),'Color','k','LineWidth',1);
%     end
%     hold on
%     plotToolbar(G.Matrix, state.pressure(1:G.Matrix.cells.num),'FaceAlpha',0.9,'EdgeAlpha',0.02)
%     title(['Reservoir Pressure - ',solvers{1,figs}],'FontSize',15,'FontWeight','bold');
%     colormap(jet)
%     colorbar
%     axis tight
% end
% %-------------------------------Surface-----------------------------------%
if isfield(G.Matrix,'cartDims')
    figure('Position',[screenSize(1:end-1),screenSize(end)/1.5])
    for figs = 1:2
        subplot(1,2,figs);
        state = getfield(sol,solvers{1,figs}); %#ok
        X = reshape(G.Matrix.cells.centroids(:,1),G.Matrix.cartDims(1),G.Matrix.cartDims(2));
        Y = reshape(G.Matrix.cells.centroids(:,2),G.Matrix.cartDims(1),G.Matrix.cartDims(2));
        P = reshape(state.pressure(1:G.Matrix.cells.num),size(X));
        surfc(X,Y,P,'EdgeAlpha',0.02,'FaceAlpha',1);
        hold on
        for i = 1:numel(fracture.lines)
            x = [fracture.lines(i).endp(1);fracture.lines(i).endp(3)];
            y = [fracture.lines(i).endp(2);fracture.lines(i).endp(4)];
            Gfl = G.FracGrid.(['Frac',num2str(i)]);
            z = [state.pressure(Gfl.cells.start);state.pressure(Gfl.cells.start+Gfl.cells.num-1)];
            line(x,y,z,'Color','k','LineWidth',1);
        end
        axis tight
        colormap(jet)
        view(30,20)
        zlabel('Pressure [Pa]','FontSize',15,'FontWeight','bold');
        title(['Reservoir Pressure - ',solvers{1,figs}],'FontSize',15,'FontWeight','bold');
        colorbar
        axis tight
    end
end
L1 = norm((state_ms.pressure(1:G.Matrix.cells.num)-...
    state_fs.pressure(1:G.Matrix.cells.num)),1)/norm(state_fs.pressure(1:G.Matrix.cells.num),1);
L2 = norm((state_ms.pressure(1:G.Matrix.cells.num)-...
    state_fs.pressure(1:G.Matrix.cells.num)),2)/norm(state_fs.pressure(1:G.Matrix.cells.num),2);
%--------------------------------Flux-------------------------------------%
% figure('Position',[screenSize(1:end-1),screenSize(end)/1.5])
% for figs = 1:2
%     subplot(1,2,figs);
%     state = getfield(sol,solvers{1,figs}); %#ok
%     cell_flux = accumarray(cellNo,...
%         abs(convertTo(faceFlux2cellFlux(G, state.flux), meter^3/day)));
%     plotToolbar(G.Matrix, cell_flux(1:G.Matrix.cells.num)); colormap(jet);
%     for i = 1:size(fl,1)
%         line(fl(i,1:2:3),fl(i,2:2:4),'Color','k','LineWidth',1);
%     end
%     title(['Flux intensity - ',solvers{1,figs}],'FontSize',15,'FontWeight','bold'); colorbar
% end

%% TOF
% figure('Position',[screenSize(1:end-1),screenSize(end)/1.5]);
% for figs = 1:2
%     subplot(1,2,figs);
%     state = getfield(sol,solvers{1,figs}); %#ok
%     hold on
%     D = computeTOFandTracer(state,G,G.rock,'Wells',W);
%     tofs = log(abs(sum(D.tof,2)));
%     m1 = max(min(tofs),mean(tofs)-2*std(tofs));
%     m2 = min(max(tofs),mean(tofs)+std(tofs));
%     plotToolbar(G.Matrix,tofs(1:G.Matrix.cells.num),'FaceAlpha',0.9,'EdgeAlpha',0.02);
%     colorbar;
%     colormap(jet);
%     caxis([m1 m2]);
%     axis tight
%     for i = 1:size(fl,1)
%         line(fl(i,1:2:3),fl(i,2:2:4),'Color','k','LineWidth',1);
%     end
%     title(['Total TOF - ',solvers{1,figs}],'FontSize',15,'FontWeight','bold'); colorbar
% end