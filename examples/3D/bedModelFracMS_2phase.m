close all;
clearvars -except METISPATH mrstVerbose screenSize GlobTri
dbstop if error

pth = getDatasetPath('bedmodel2');

grdecl = readGRDECL(fullfile(pth, 'BedModel2.grdecl'));
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

G = processGRDECL(grdecl);
G = computeGeometry(G);
% grdecl.PERMX(grdecl.SATNUM==3) = 1e-20;
G.rock = grdecl2Rock(grdecl, G.cells.indexMap);
load('GlobTri_BedModel2.mat');
% global GlobTri
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
fluid = initSimpleFluid('mu' , [   1,  1] .* centi*poise     , ...
    'rho', [1000, 700] .* kilogram/meter^3, ...
    'n'  , [   2,   2]);

%%
M = G;
G = assembleGlobalGrid(G); 
pv = abs(poreVolume(G,G.rock));
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
pinit = 0;
state = initResSol (G, pinit, 0);
% Get A matrix without source
A = incompTPFA(state, G, T, fluid, 'MatrixOutput', true, 'use_trans',true); 
A = A.A; A(1,1) = A(1,1)/2; % undo magic done in incompTPFA

%% Setup Multiscale Grids
G.fractureAperture = fracplanes(1).aperture;
dispif(mrstVerbose, 'RSB Grids...\n\n');
dof_frac = [3 3 2 2 2];
nw = struct; 
for i = 1:numel(fracplanes)
    nw(i).lines = i;
end
% 
CG = getRsbGrids_HFM(G, nw, 'coarseDims', [6 6 5], 'sampleDims', [30 30 60], 'dof_frac',dof_frac, 'sysMatrix', A);
dof_matrix = max(CG.partition(1:G.Matrix.cells.num));

%% Interaction Regions
IR = NaN(G.cells.num,CG.cells.num);
for i = 1:CG.cells.num
IR(CG.cells.interaction{i,1},i) = 1;
end
figure; plotToolbar(G, IR);
figure; plotToolbar(G, IR, G.Matrix.cells.num+1:G.cells.num);

%%
W = [];
xf = G.faces.centroids(:, 1);
left = find(abs(xf - min(xf)) < 1e-4);
right = find(abs(xf - max(xf)) < 1e-4);

bc = [];
bc = addBC(bc, left, 'pressure', 1, 'sat', [1 0]);
bc = addBC(bc, right, 'pressure', 0, 'sat', [0 1]);


%% Init state

state_fs = incompTPFA(state, G, T, fluid,  ...
    'bc', bc, 'Wells', W, 'MatrixOutput', true, 'use_trans',true);

% Multiscale

basis_sb = getMultiscaleBasis(CG, A, 'type', 'rsb');
state_ms = incompMultiscale(state, CG, T, fluid, basis_sb, 'Wells', W, ...
    'bc', bc,'use_trans',true);

%% Transport
pv     = poreVolume(G,G.rock);
nt     = 15;
t90    = 0.9*(sum(pv)/abs(sum(state_fs.flux(left))));
Time   = t90;
dT     = Time/nt;

sol_fs = cell(nt,1); sol_ms = cell(nt,1); Bfrac = cell(nt,1); Bfull = Bfrac;
pvi = zeros(nt,1);
timesteps = repmat(dT, nt, 1); timesteps = cumsum(timesteps);
L1 = zeros(nt,1); L2 = zeros(nt,1); L_inf = zeros(nt,1);
% state_ms.rhs = rhs;
t  = 0; plotNo = 1; e = []; 
dispif(mrstVerbose, '\nTransport...\n\n');
count = 1;
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
    B = iteratedJacobiBasis(A, CG, 'interpolator', basis_sb.B); 
    R = controlVolumeRestriction(CG.partition);
    basis_sb = struct('B', B, 'R', R);
    basis_frac = struct('B',B(G.Matrix.cells.num+1:end,:));
    Bfrac{count,1} = basis_frac;
    Bfull{count,1} = basis_sb;
    state_ms = incompMultiscale(state_ms, CG, T, fluid, basis_sb, 'Wells', W,'use_trans',true);
    sol_fs{count,1} = state_fs; sol_ms{count,1} = state_ms; 
    
    t = t + dT;
    pvi(count) = 100*(sum(state_fs.flux(left))*t)/sum(pv);
%     [state_ms_iter,report] = incompMultiscale(state_ms, CG, T, fluid, basis_sb,...
%     'Wells', W, 'bc', bc, 'use_trans',true, 'tolerance', 1e-6, 'iterations', 1e2,...
%     'useGMRES', false, 'reconstruct', true, 'getSmoother', fn);

    L1(count) = norm((state_ms.s-state_fs.s).*pv,1)/norm(state_fs.s.*pv,1);
    L_inf(count) = norm((state_ms.s-state_fs.s).*pv,'inf')/norm(state_fs.s.*pv,'inf');
    L2(count) = sum(((state_ms.s-state_fs.s).^2) .* pv)/sum((state_fs.s.^2).*pv);
    % Increase time and continue if we do not want to plot saturations
    fprintf([num2str(pvi(count)), '%% PVI \n']);
    count = count+1;
end
   
%%
    
figure; plotToolbar(G,sol_fs);
figure; plotToolbar(G,sol_ms);
figure; plotToolbar(CGf.parent,basis_frac)
figure; plot(pvi,L1)
eq = ['$$ \mathbf{ \frac{\sum_{i} | S_i^{fs}-S_i^{ms} | \times pv_i}{\sum_{i} pv_i}  \times 100 } $$'];
title([eq, '  \textbf{vs PVI}'],...
        'FontSize',16,'FontWeight', 'bold','interpreter','latex');
ax = gca;
ax.YLabel.String = 'L_1';
ax.XLabel.String = '% PVI';