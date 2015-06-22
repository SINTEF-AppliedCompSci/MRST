close all;
clearvars -except METISPATH mrstVerbose screenSize
% dbstop if warning
dbstop if error

% global GlobTri
nx = 20; ny = 20; nz = 20;
lx = 100; ly = 100; lz = 100;
G = cartGrid([nx ny nz], [lx ly lz]);
G = computeGeometry(G);
% if isempty(GlobTri)
%     GlobTri = globalTriangulation(G,GlobTri);
% end

%%
fracplanes = planes_input();
checkIfCoplanar(fracplanes)

%%
[G,fracplanes] = preProcessingFractures(G, fracplanes);
fprintf('\nPreprocessing complete.\n');

%%
G.rock.perm = ones(G.cells.num,1)*darcy();
G.rock.poro = 0.2*ones(G.cells.num, 1);
K_frac = 10000; % Scaling factor = K_frac/K_mat D
poro_frac = 0.5;
G = makeRockFrac(G, K_frac, 'permtype','homogeneous','rockporo',poro_frac);
fluid = initSingleFluid('mu' , 1*centi*poise, ...
    'rho', 1000*kilogram/meter^3);

%%
% G.fractureAperture = fracture.aperture;
G = assembleGlobalGrid(G);
w1 = G.cells.volumes(G.nnc.cells(:,1))./G.rock.perm(G.nnc.cells(:,1));
w2 = G.cells.volumes(G.nnc.cells(:,2))./G.rock.perm(G.nnc.cells(:,2));
wt = G.cells.volumes(G.nnc.cells(:,1))+G.cells.volumes(G.nnc.cells(:,2));
G.nnc.T = G.nnc.CI.*(wt./(w1+w2)); clear wt w1 w2


%% compute transmissibilities
T = computeTrans(G, G.rock);
%-------------------------------------------------------------------------%
% computeTrans returns 2 transmissibilities for each internal face and one
% transmissibility fo each external face. Size of transmissibility array is
% same as G.cells.faces if opt.usetrans = false 
%-------------------------------------------------------------------------%
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
T = [T;G.nnc.T];

%%
pinit = 100*barsa();
state = initResSol (G, pinit);
wpinit = 100*barsa();
state.wellSol = initWellSol(G, wpinit);
% Get A matrix without source
A = incompTPFA(state, G, T, fluid, 'MatrixOutput', true, 'use_trans',true); 
A = A.A; A(1,1) = A(1,1)/2; % undo magic done in incompTPFA

%% Setup Multiscale Grids
G.fractureAperture = 1/25;
dispif(mrstVerbose, 'RSB Grids...\n\n');
coarsen = [5 5 5]; % coarsening factor in each direction
dof_frac = [2 2];
dof_matrix = 100; % Coarse dof
nw = struct; 
for i = 1:numel(fracplanes)
    nw(i).lines = i;
end
% set partition_frac to false if you don't want a separate interaction
% region for fractures. i.e. no rsb in fractures
[CG, CGf] = getRsbGrids_HFM(G, nw, 'coarsen', coarsen, ...
    'use_metis', false,'dof_matrix',dof_matrix,'dof_frac',dof_frac,'sysMatrix',A);
dof_matrix = max(CG.partition(1:G.Matrix.cells.num));


%%
radius = 0.1;
% cellinj = 1:nx*ny:nx*ny*nz;
% cellprod = nx*ny:nx*ny:nx*ny*nz;
cellinj = nx*ny*(nz-1)+1:nx:nx*ny*nz;
cellprod = nx:nx:nx*ny;
W   = addWell([], G.Matrix, G.Matrix.rock, cellinj, 'Type', 'rate',...
    'Val', 1000/day(), 'Radius', radius, 'Sign',1);
W   = addWell(W, G.Matrix, G.Matrix.rock, cellprod, ...
    'Type', 'bhp', 'Val', 100*barsa(), 'Radius', radius, 'Sign', -1);


%% Incompressible 1-phase FS

bc = []; %
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
figure; plotToolbar(G,basis_sb.B); colormap(jet);
%----------------------------1 MSFV operation-----------------------------%

[state_ms,report1] = incompMultiscale(state, CG, T, fluid, basis_sb, 'Wells', W, ...
    'bc', bc,'use_trans',true);

%----------------------------Use ILU Smoother-----------------------------%
fn = getSmootherFunction('type', 'ilu');

[state_ms_iter,report] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'Wells', W, 'bc', bc, 'use_trans',true, 'tolerance', 1e-6, 'iterations', 1e3,...
    'useGMRES', false, 'reconstruct', true, 'getSmoother', fn);
sol = {state_fs;state_ms}; 


%% Plots - pressure

dispif(mrstVerbose, 'Plots...\n\n');
figure
%-------------------------------Patch-------------------------------------%
plotToolbar(G, sol,'FaceAlpha',0.9,'EdgeAlpha',0.02)
plotWell(G,W); view(-20,20);
colormap(jet)
colorbar
axis tight

L1 = norm((state_ms.pressure(1:G.Matrix.cells.num)-...
    state_fs.pressure(1:G.Matrix.cells.num)),1)/norm(state_fs.pressure(1:G.Matrix.cells.num),1);
L_inf = norm((state_ms.pressure-...
    state_fs.pressure),'inf')/norm(state_fs.pressure,'inf');
pv     = poreVolume(G,G.rock);
L2 = sum(((state_ms.pressure-...
    state_fs.pressure).^2) .* pv)/sum((state_fs.pressure.^2).*pv);


ploterr = abs(state_ms.pressure-...
    state_fs.pressure)./state_fs.pressure;
figure;
plotToolbar(G, ploterr,'FaceAlpha',0.9,'EdgeAlpha',0.02)
plotWell(G,W); view(-20,20);
locale_eq = '$$ \mathbf{ \frac{| P_i^{fs}-P_i^{ms} | }{ | P_i^{fs} |} } $$';
title(locale_eq,'FontSize',15,'FontWeight','bold','interpreter','latex');
colormap(jet)
colorbar
axis tight