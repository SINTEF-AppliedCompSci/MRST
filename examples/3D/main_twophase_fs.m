close all;
clearvars -except METISPATH mrstVerbose screenSize
% dof_frac = 30; % Coarse dof
% dof_matrix = 100; % Coarse dof
% dbstop if warning
dbstop if error

% global GlobTri
nx = 20; ny = 20; nz = 20;
lx = 100; ly = 100; lz = 100;
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
fluid = initSimpleFluid('mu' , [   1,  1] .* centi*poise     , ...
    'rho', [1000, 700] .* kilogram/meter^3, ...
    'n'  , [   2,   2]);

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
pinit = 175*barsa();
state = initResSol (G, pinit, 0);
wpinit = 150*barsa();
state.wellSol = initWellSol(G, wpinit);

radius = 0.1;
% cellinj = 1:nx*ny:nx*ny*nz;
% cellprod = nx*ny:nx*ny:nx*ny*nz;
cellinj = nx*ny*(nz-1)+1:nx:nx*ny*nz;
cellprod = nx:nx:nx*ny;
W   = addWell([], G.Matrix, G.Matrix.rock, cellinj, 'Type', 'bhp',...
    'Val', 200*barsa(), 'Radius', radius, 'Sign',1, 'Comp_i', [1, 0]);
W   = addWell(W, G.Matrix, G.Matrix.rock, cellprod, ...
    'Type', 'bhp', 'Val', 150*barsa(), 'Radius', radius, 'Sign', -1, 'Comp_i', [0, 1]);

bc = [];
dispif(mrstVerbose, 'Computing initial reservoir state...\n\n');
state = incompTPFA(state, G, T, fluid,  ...
    'bc',bc, 'Wells', W, 'MatrixOutput', true, 'use_trans',true);


%% Transport
pv     = poreVolume(G,G.rock);
nt     = 15;
t90   = 0.9*(sum(pv)/abs(sum(state.wellSol(1).flux)));
Time   = t90;
dT     = Time/nt;

pvi = zeros(nt,1);
sol = cell(nt,1);
t  = 0; 
dispif(mrstVerbose, '\nTransport...\n\n');
count = 1;
while t < Time,
    state = implicitTransport(state, G, dT, G.rock, fluid, 'wells', W, 'bc', bc, 'Trans', T);
    % Check for inconsistent saturations
    s = state.s(:,1);
    assert(max(s) < 1+eps && min(s) > -eps);

    % Update solution of pressure equation.
    state  = incompTPFA(state , G, T, fluid, 'wells', W, 'bc', bc, 'use_trans',true);
    
    t = t + dT;
    pvi(count) = 100*(sum(state.wellSol(1).flux)*t)/sum(pv);
    fprintf([num2str(pvi(count)), '%% PVI \n']);
    sol{count,1} = state;
    count = count+1;
end
   
%%
wellSol = convertIncompWellSols(W, sol);
    
figure; plotToolbar(G,sol); plotWell(G,W); view(-20,20);
plotWellSols(wellSol,pvi); 
%%
[neighborship, nnc] = getNeighbourship(G,'Topological',1);
[cellNo, cellFaces, isNNC] = getCellNoFaces(G);
nncfaces = cellFaces(isNNC);
flux_nnc = sum(state_fs.flux(nncfaces));
cell_flux = accumarray(cellNo,...
        abs(convertTo(faceFlux2cellFlux(G, state.flux), meter^3/day)));
    plotToolbar(G.Matrix, cell_flux(1:G.Matrix.cells.num)); colormap(jet);

%%
fn = @(x) sqrt(sum(faceFlux2cellVelocity(G, x.flux).^2, 2));
v = cellfun(fn, sol, 'uniformoutput', false);
figure; plotToolbar(G, v)