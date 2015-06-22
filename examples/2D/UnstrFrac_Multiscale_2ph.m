% close all;
clearvars -except METISPATH mrstVerbose screenSize
fgridtype = 1; % approx. cell size specified
min_size = 2; cell_size = 5; plotperm = 0;
a = 1/25; % fracture aperture
dof_frac = 5; % Coarse dof
dof_matrix = 49; % Coarse dof
fracIRtype = 4; % connectivity based
%% User defined fracture lines and appropriate grid
load('Frac_lines_10_5nw.mat');
load('Gtri_100x100.mat');
% G = pebi(G);
G = computeGeometry(G);


%% Process fracture lines into grid

[G,fracture] = processFracture2D(G,fl); fracture.aperture = a;
figure; plotFractureLines(G,fracture,'network'); % 3rd arg can be line #'s, network or none
% figure; plotMarkedCells(G, fracture,'lines'); % 3rd arg can be line #'s, network or none

%% Compute CI and grid fracture

dispif(mrstVerbose, 'Computing CI and gridding fracture...\n\n');
G = CIcalculator2D(G,fracture);
% plotCI(G,fracture);

[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
% figure; plotFractureNodes2D(G,F,fracture,show_node_num); clear show_node_num
%% Set rock permeabilities in fracture and matrix

dispif(mrstVerbose, 'Initializing rock properties...\n\n');
G.rock.perm = ones(G.cells.num,1)*darcy();
% G.rock.perm = randPermGen(G.cells.num);
% G.rock.perm = randi(5e10,G.cells.num,1)*darcy()./randi(5e12,G.cells.num,1);
G.rock.poro = 0.2*ones(G.cells.num, 1);
K_frac = 10000; % Scaling factor = K_frac/K_mat D
poro_frac = 0.5;
G = makeRockFrac(G, K_frac, 'permtype','homogeneous','rockporo',poro_frac);
if plotperm == 1, plotMatrixPermeability(G);  end


%% Define NNC and corresponding Trans

[G,T] = defineNNCandTrans(G,F,fracture);

%% Define the two-phase fluid model
%
% * densities: [rho_w, rho_o] = [1000 700] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 10] cP.
fluid = initSimpleFluid('mu' , [   1,  1] .* centi*poise     , ...
    'rho', [1000, 700] .* kilogram/meter^3, ...
    'n'  , [   1,   1]);

%% Relperms

s = linspace(0,1,101)'; kr=fluid.relperm(s);
% figure;
% plot(s, kr(:,1), 'b', s, kr(:,2), 'r');
% title('Relative permeability curves')
% legend('Water','Oil','Location','Best')


%% Initialize state variables

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
pinit = 5*barsa(); wpinit = 5*barsa();
state  = initResSol (G, pinit);
state.wellSol = initWellSol(G, wpinit);
% Get A matrix without source
A = incompTPFA(state, G, T, fluid, 'MatrixOutput', true, 'use_trans',true); 
rhs = A.rhs; A = A.A; A(1,1) = A(1,1)/2; % undo magic done in incompTPFA


%% Setup Multiscale Grids

dispif(mrstVerbose, 'Defining coarse grids and interaction regions...\n\n');

coarsen = [8 8]; % coarsening factor in each direction
 
% set partition_frac to false if you don't want a separate interaction
% region for fractures. i.e. no rsb in fractures
[CG, CGf] = getRsbGrids_HFM(G, fracture.network, 'coarsen', coarsen, ...
    'use_metis', true, 'dof_matrix', dof_matrix, 'dof_frac', dof_frac,...
    'sysMatrix', A, 'fracIRtype', fracIRtype, 'fullyCoupled', true);

dof_matrix = max(CG.partition(1:G.Matrix.cells.num));
figure; plotFractureCoarseGrid2D(G,CG.partition,F)


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
% figure; plotFracData2D(G, IRf, 'figTitle','Internal fracture interaction regions',...
%     'wide',true,'F',F,'CG',CGf,'plotMatrix',true,'cmap',flipud(winter));


%% Add Wells and BC

% No bc = no-flow boundaries
bc = [];
wellRadius = 9.5e-4;
wcells = findQuarter5spotLocations(G,[0 0],[100 100]);
W = addWell([], G.Matrix, G.Matrix.rock, wcells(1),'InnerProduct', 'ip_tpf','Type', ...
    'rate', 'Val', 100/day(),'Radius', wellRadius, 'Comp_i', [1, 0]);
W = addWell(W, G.Matrix, G.Matrix.rock, wcells(2), 'InnerProduct', 'ip_tpf', 'Type', ...
    'bhp' , 'Val', 10*barsa(), 'Radius', wellRadius, 'Comp_i', [0, 1]);

% Generate the components of the linear system corresponding to the two
% wells (or bc) and initialize the solution structure (with correct bhp)
% initState can only be called when wells are used
statei = initState(G, W, 0, [0, 1]); % Initialize State

%% Solve for initial pressure in the reservoir

dispif(mrstVerbose, 'Computing initial reservoir state...\n\n');
state_fs = incompTPFA(state, G, T, fluid,  ...
    'bc',bc, 'Wells', W, 'MatrixOutput', true, 'use_trans',true);
[neighborship, nnc] = getNeighbourship(G,'Topological',1);
[cellNo, cellFaces, isNNC] = getCellNoFaces(G);
nncfaces = cellFaces(isNNC);
flux_nnc = sum(state_fs.flux(nncfaces));


%% Multiscale basis and solution

basis_sb = getMultiscaleBasis(CG, A, 'type', 'rsb', 'useMex', false);
% figure; plotToolbar(G,basis_sb.B); axis tight; c = colormap(jet);  
% c(1,:) = [1 1 1]; colormap(c); colorbar; 
% title('Basis Functions','FontSize',15,'FontWeight','bold');

%-----------Dirty way to see basis functions inside fractures-------------%

% figure; plotFracData(G, basis_sb.B(G.Matrix.cells.num+1:G.cells.num,dof_matrix+1:max(CG.partition)), ...
%     'figTitle','Basis function inside fractures',...
%     'wide',true,'F',F,'CG',CGf,'plotMatrix',true,'cmap',jet);

%----------------------------1 MSFV operation-----------------------------%

[state_ms,report1] = incompMultiscale(statei, CG, T, fluid, basis_sb, 'Wells', W, ...
    'bc', bc,'use_trans',true);

%----------------------------Use ILU Smoother-----------------------------%
fn = getSmootherFunction('type', 'ilu');

[state_ms_iter,report] = incompMultiscale(statei, CG, T, fluid, basis_sb,...
    'Wells', W, 'bc', bc, 'use_trans',true, 'tolerance', 1e-6, 'iterations', 1e3,...
    'useGMRES', false, 'reconstruct', true, 'getSmoother', fn);
sol = struct('FS',state_fs,'MS',state_ms); solvers = {'FS', 'MS'};


%% Report initial reservoir pressure

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
%-------------------------------Surface-----------------------------------%
% if isfield(G.Matrix,'cartDims')
%     figure('Position',[screenSize(1:end-1),screenSize(end)/1.5])
%     for figs = 1:2
%         subplot(1,2,figs);
%         state = getfield(sol,solvers{1,figs}); %#ok
%         X = reshape(G.Matrix.cells.centroids(:,1),G.Matrix.cartDims(1),G.Matrix.cartDims(2));
%         Y = reshape(G.Matrix.cells.centroids(:,2),G.Matrix.cartDims(1),G.Matrix.cartDims(2));
%         P = reshape(state.pressure(1:G.Matrix.cells.num),size(X));
%         surfc(X,Y,P,'EdgeAlpha',0.02,'FaceAlpha',1);
%         hold on
%         for i = 1:numel(fracture.lines)
%             x = [fracture.lines(i).endp(1);fracture.lines(i).endp(3)];
%             y = [fracture.lines(i).endp(2);fracture.lines(i).endp(4)];
%             Gfl = G.FracGrid.(['Frac',num2str(i)]);
%             z = [state.pressure(Gfl.cells.start);state.pressure(Gfl.cells.start+Gfl.cells.num-1)];
%             line(x,y,z,'Color','k','LineWidth',1);
%         end
%         axis tight
%         colormap(jet)
%         view(30,20)
%         zlabel('Pressure [Pa]','FontSize',15,'FontWeight','bold');
%         title(['Initial Reservoir Pressure - ',solvers{1,figs}],'FontSize',15,'FontWeight','bold');
%         colorbar
%         axis tight
%     end
% end
% %--------------------------------Flux-------------------------------------%
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
%     title(['Initial flux intensity - ',solvers{1,figs}],'FontSize',15,'FontWeight','bold'); colorbar
% end

%% Incompressible Two-Phase Flow
%
%---------------------------Transport loop--------------------------------%
% We solve the two-phase system using a sequential splitting in which the
% pressure and fluxes are computed by solving the flow equation and then
% held fixed as the saturation is advanced according to the transport
% equation. This procedure is repeated for a given number of time steps The
% error introduced by this splitting of flow and transport can be reduced
% by iterating each time step until e.g., the residual is below a certain
% user-prescribed threshold (this is not done herein).
pv     = poreVolume(G,G.rock);
nt     = 30;
t90    = 0.9*(sum(pv)/state_fs.wellSol(1).flux);
Time   = t90;
dT     = Time/nt;
dTplot = Time/3;  % plot only every 1500 days
N      = fix(Time/dTplot);

sat_error = zeros(G.Matrix.cells.num,Time/dTplot); L1 = zeros(Time/dTplot,1);
state_ms.rhs = rhs;
t  = 0; plotNo = 1; hfs = 'Fine-scale: '; hms = 'Multiscale: ';
e = []; pms = []; pfs = [];
dispif(mrstVerbose, 'Solving Transport...\n\n');
figure
B = basis_sb.B;
while t < Time,
    state_fs = implicitTransport(state_fs, G, dT, G.rock, fluid, 'wells', W, 'bc', bc, 'Trans', T);
    state_ms = implicitTransport(state_ms, G, dT, G.rock, fluid, 'wells', W, 'bc', bc, 'Trans', T);
    % Check for inconsistent saturations
    s = [state_fs.s(:,1); state_ms.s(:,1)];
    assert(max(s) < 1+eps && min(s) > -eps);

    % Update solution of pressure equation.
    state_fs  = incompTPFA(state_fs , G, T, fluid, 'wells', W, 'bc', bc, 'use_trans',true);

    %-------------------------------Multiscale----------------------------%
    A = incompTPFA(state_ms, G, T, fluid, 'MatrixOutput', true, ...
        'use_trans',true); A = A.A;
    B = iteratedJacobiBasis(A, CG, 'interpolator', B); 
    R = controlVolumeRestriction(CG.partition);
    basis_sb = struct('B', B, 'R', R);
    state_ms = incompMultiscale(state_ms, CG, T, fluid, basis_sb, 'Wells', W,'use_trans',true);

    [state_ms_iter,report] = incompMultiscale(state_ms, CG, T, fluid, basis_sb,...
    'Wells', W, 'bc', bc, 'use_trans',true, 'tolerance', 1e-6, 'iterations', 1e2,...
    'useGMRES', false, 'reconstruct', true, 'getSmoother', fn);

    %---------------------------------------------------------------------%

    % Measure water saturation in production cells in saturation
    e = [e; sum(abs(state_fs.s(:,1) - state_ms.s(:,1)).*pv)/sum(pv)]; %#ok
    pfs = [pfs; state_fs.s(W(2).cells,1)' ];                 %#ok
    pms = [pms; state_ms.s(W(2).cells,1)'];                 %#ok
    
    % Increase time and continue if we do not want to plot saturations
    t = t + dT;
    if ( t < plotNo*dTplot && t <Time), continue, end
    
    % Plot saturation
    pvi = 100*(state_fs.wellSol(1).flux*t)/sum(pv);
    heading = [num2str(round(pvi,1)),  ' % PVI'];
    r = 0.01;
    subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.48]), cla
    plotToolbar(G.Matrix, state_fs.s(1:G.Matrix.cells.num,1),'FaceAlpha',0.9);
    for i = 1:size(fl,1)
        line(fl(i,1:2:3),fl(i,2:2:4),'Color','k','LineWidth',1);
    end
    axis equal off, title([hfs heading],'FontSize',14,'FontWeight','bold');
    colormap(flipud(jet)), view(0,90)
    
    subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.48]), cla
    plotToolbar(G.Matrix, state_ms.s(1:G.Matrix.cells.num,1),'FaceAlpha',0.9);
    for i = 1:size(fl,1)
        line(fl(i,1:2:3),fl(i,2:2:4),'Color','k','LineWidth',1);
    end
    axis equal off, title([hms heading],'FontSize',14,'FontWeight','bold');
    colormap(flipud(jet)), view(0,90)
    sat_error(:,plotNo) = abs(state_ms.s(1:G.Matrix.cells.num,1)...
        - state_fs.s(1:G.Matrix.cells.num,1));
    L1(:,plotNo) = sum(abs(state_ms.s(1:G.Matrix.cells.num,1)...
        - state_fs.s(1:G.Matrix.cells.num,1)))/sum(state_fs.s(1:G.Matrix.cells.num,1));
    plotNo = plotNo+1;
end
   
%% Plot pore volume weighted L1 norm for saturation and error at producer
%
figure('Position',[screenSize(1:end-1),screenSize(end)/1.5])
n = size(pfs,1);
subplot(1,2,1),
l1_eq = '$$ \mathbf{ \frac{\sum_{i} | S_i^{fs}-S_i^{ms} | \times pv_i}{\sum_{i} pv_i}  \times 100}$$';
plot(1:n,e*100,'-o'), title([l1_eq, '  \textbf{vs Timestep}'],...
    'FontSize',16,'FontWeight', 'bold','interpreter','latex');
ylabel('%','FontSize',16,'FontWeight', 'bold');
xlabel('Time Step','FontSize',16,'FontWeight', 'bold'); axis tight
subplot(1,2,2),
plot(1:n,pfs(:,1),'-o',1:n,pms(:,1),'--*')
leg = legend('Fine-scale','Multiscale','Location','Best');
set(leg,'FontSize',16,'FontWeight', 'bold');
title('Water breakthrough (saturation) at producer', 'FontSize',16,'FontWeight', 'bold');
ylabel('Saturation at producer','FontSize',16,'FontWeight', 'bold');
xlabel('Time Step','FontSize',16,'FontWeight', 'bold'); axis tight