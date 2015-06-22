close all;
clearvars -except METISPATH mrstVerbose screenSize
% dof_frac = 30; % Coarse dof
% dof_matrix = 100; % Coarse dof
% dbstop if warning
dbstop if error

% global GlobTri
nx = 20; ny = 20; nz = 5;
lx = 100; ly = 100; lz = 20;
G = cartGrid([nx ny nz], [lx ly lz]);
G = computeGeometry(G);
G.Domain = [lx ly lz];
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
% wpinit = 100*barsa();
% state.wellSol = initWellSol(G, wpinit);

radius = 0.1;
% cellinj = 1:nx*ny:nx*ny*nz;
% cellprod = nx*ny:nx*ny:nx*ny*nz;
cellinj = nx*ny*(nz-1)+1:nx:nx*ny*nz;
cellprod = nx:nx:nx*ny;
W   = addWell([], G.Matrix, G.Matrix.rock, cellinj, 'Type', 'rate',...
    'Val', 1000/day(), 'Radius', radius, 'Sign',1);
W   = addWell(W, G.Matrix, G.Matrix.rock, cellprod, ...
    'Type', 'bhp', 'Val', 100*barsa(), 'Radius', radius, 'Sign', -1);

bc = [];
dispif(mrstVerbose, 'Computing initial reservoir state...\n\n');
state = incompTPFA(state, G, T, fluid,  ...
    'bc',bc, 'Wells', W, 'MatrixOutput', true, 'use_trans',true);
figure; plotToolbar(G,state); plotWell(G,W);
caxis([min(state.pressure) max(state.pressure)]);