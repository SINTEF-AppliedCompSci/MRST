mrstModule add new-multiscale coarsegrid ad-core
close all;
clearvars -except METISPATH mrstVerbose screenSize GlobTri
% dbstop if error

pth = getDatasetPath('bedmodel2');

grdecl = readGRDECL(fullfile(pth, 'BedModel2.grdecl'));
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

G = processGRDECL(grdecl);
G = computeGeometry(G);

G.rock = grdecl2Rock(grdecl, G.cells.indexMap);
load('GlobTri_BedModel2.mat');
global GlobTri
% if isempty(GlobTri)
%     GlobTri = globalTriangulation(G,GlobTri);
% elseif size(GlobTri.Tri.Points,1) ~= G.nodes.num
%     GlobTri = globalTriangulation(G,GlobTri);
% end

%%
fracplanes = planes_input();
checkIfCoplanar(fracplanes)

%%
[G,fracplanes] = preProcessingFractures(G, GlobTri, fracplanes);
fprintf('\nPreprocessing complete.\n');

%%
% G.rock.perm = ones(G.cells.num,1);
K_frac = 10000; % Scaling factor = K_frac/K_mat D
poro_frac = 0.5;
G = makeRockFrac(G, K_frac, 'permtype','homogeneous','rockporo',poro_frac);
fluid = initSingleFluid('mu' , 1*centi*poise, ...
    'rho', 1000);

%%
% G.fractureAperture = fracture.aperture;
M = G;
G = assembleGlobalGrid(G); 
pv = poreVolume(G,G.rock);
w1 = pv(G.nnc.cells(:,1))./G.rock.perm(G.nnc.cells(:,1));
w2 = pv(G.nnc.cells(:,2))./G.rock.perm(G.nnc.cells(:,2));
wt = pv(G.nnc.cells(:,1)) + pv(G.nnc.cells(:,2));
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
state = initResSol(G, 0);
A = getIncomp1PhMatrix(G, T, state, fluid);

%% Setup Multiscale Grids
G.fractureAperture = 1/25;
dispif(mrstVerbose, 'RSB Grids...\n\n');
coarsen = [6 6 10]; % coarsening factor in each direction
dof_frac = 4;
dof_matrix = 100; % Coarse dof
nw = struct; 
for i = 1:numel(fracplanes)
    nw(i).lines = i;
end

%% Grid of roughly same dimensions
dims = max(G.Matrix.nodes.coords);
Gc = cartGrid([30, 30, 60], dims);
Gc = computeGeometry(Gc);
%% Make coarse grid
p = partitionUI(Gc, [6, 6, 5]);
p = reshape(p, Gc.cartDims);
p = sampleFromBox(G.Matrix, p);
%%
CG = olavFractureCoarsen(G, p, 4, 5);

%%
W = [];
xf = G.faces.centroids(:, 1);
left = find(abs(xf - min(xf)) < 1e-4);
right = find(abs(xf - max(xf)) < 1e-4);

bc = [];
bc = addBC(bc, left, 'pressure', 1);
bc = addBC(bc, right, 'pressure', 0);


%% Incompressible 1-phase FS

dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state_ref = incompTPFA(state, G, T, fluid,  ...
    'bc', bc , 'Wells', W, 'MatrixOutput', true, 'use_trans',true);

%% Multiscale basis and solution
mrstModule add mrst-gui

dispif(mrstVerbose, 'Solving multiscale system...\n\n');

%----------------------------Basis functions------------------------------%

basis_sb = getMultiscaleBasis(CG, A, 'type', 'rsb', 'useControlVolume' , true, 'useMex', true);
figure; plotToolbar(G,basis_sb.B); colormap(jet);
%----------------------------1 MSFV operation-----------------------------%

[state_ms,report1] = incompMultiscale(state, CG, T, fluid, basis_sb, 'Wells', W, ...
    'bc', bc,'use_trans',true);

%----------------------------Use ILU Smoother-----------------------------%
fn = getSmootherFunction('type', 'ilu');

[state_ms_iter,report] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'Wells', W, 'bc', bc, 'use_trans',true, 'tolerance', 1e-6, 'iterations', 1e3,...
    'useGMRES', false, 'reconstruct', true, 'getSmoother', fn);
sol = {state_ref;state_ms}; 


%% Plots - pressure

dispif(mrstVerbose, 'Plots...\n\n');
figure
%-------------------------------Patch-------------------------------------%
plotToolbar(G, state_ref,'FaceAlpha',0.9,'EdgeAlpha',0.02)
colormap(jet)
colorbar
axis tight

figure
plotToolbar(G, state_ms,'FaceAlpha',0.9,'EdgeAlpha',0.02)
caxis([0, 1])
colormap(jet)
colorbar
axis tight

L1 = norm((state_ms.pressure-state_ref.pressure),1)/norm(state_ref.pressure,1);
L_inf = norm((state_ms.pressure-state_ref.pressure),'inf')/norm(state_ref.pressure,'inf');
L2 = sum(((state_ms.pressure-state_ref.pressure).^2) .* pv)/sum((state_ref.pressure.^2).*pv);


ploterr = abs(state_ms.pressure-state_ref.pressure);
figure;
plotToolbar(G, ploterr,'FaceAlpha',0.9,'EdgeAlpha',0.02)
locale_eq = '$$ \mathbf{ \frac{| P_i^{fs}-P_i^{ms} | }{ | P_i^{fs} |} } $$';
title(locale_eq,'FontSize',15,'FontWeight','bold','interpreter','latex');
colormap(jet)
colorbar
axis tight
