%% Water injection into a 3D fractured porous media
% Two-phase example with a horizontal producer and injector simulating
% water injection in a 3-dimensional fractured porous media using the HFM
% module. Note that the 3D solvers are not capable of handling intersecting
% fracture planes.

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add coarsegrid;      % functionality for coarse grids
mrstModule add ad-core;         % NNC support for coarse grids
mrstModule add msrsb;           % MsRSB solvers
mrstModule add mrst-gui;        % plotting routines
mrstModule add incomp;          % Incompressible fluid models
checkMATLABversionHFM;

%% Grid and fracture(s)
% Construct a Cartesian grid comprising 25-by-25-by-25 cells, where each
% cell has dimension 4-by-4-by-4 m^3. Define a fracture plane using a
% coplanar set of points supplied as rows in the matrix [X Y Z]. In this
% example, we define a fracture plane of size 50-by-50 m^2, extending in
% the y-direction through the centre of the domain. Fracture aperture is
% set to 0.02 meters. Additionally, we compute a triangulation using the
% matrix grid nodes. The triangulation is used to locate cells penetrated
% by the fracture plane in the next section.
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
% Using the input fracture, we identify fine-cells in the matrix containing
% these fractures. Next, a global conductivity index (CI) is computed for
% the fracture plane and a grid is defined over it. A 'non-neighboring
% connection (NNC)' status is given to each fracture-matrix connection and
% we compute the transmissibility for each NNC using the global CI.
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
% Set the permeability (K) as 1 Darcy in the matrix and 10000 Darcy in the
% fractures. Additionally, set the porosity of the matrix and fractures to
% 20% and 50%, respectively.
dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
G.rock.perm = ones(G.cells.num,1)*darcy();
G.rock.poro = 0.2*ones(G.cells.num, 1);
K_frac = 10000; % Darcy
poro_frac = 0.5;
G = makeRockFrac(G, K_frac, 'permtype','homogeneous','porosity', poro_frac);

%% Define fluid properties
% Define a two-phase fluid model without capillarity. The fluid model has
% the values:
%
% * densities: [rho_w, rho_o] = [1000 700] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 1] cP.
% * corey-coefficient: [2, 2] = [2 2].
fluid = initSimpleFluid('mu' , [   1,  1] .* centi*poise     , ...
    'rho', [1000, 700] .* kilogram/meter^3, ...
    'n'  , [   2,   2]);

%% Assemble global grid and compute transmissibilities
% In this section, we combine the fracture and matrix grids into one grid.
% The transmissibility for each face in the combined grid and each NNC is
% computed and stored in the vector T.
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

%% Add wells
% We will include two horizontal wells, a rate-controlled injector well
% through the bottom left corner of the grid and a producer, controlled by
% bottom-hole pressure, through the top right section of the grid. Wells
% are described using a Peaceman model, giving an extra set of equations
% that need to be assembled, see simpleWellExample.m for more details on
% the Peaceman well model.
[nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
cellinj = nx*ny*(nz-1)+1:nx:nx*ny*nz;
cellprod = nx:nx:nx*ny;
W   = addWell([], G.Matrix, G.Matrix.rock, flipud(cellinj), 'Type', 'rate',...
    'Val', 500/day, 'Sign',1, 'Comp_i', [1, 0], 'Name', 'Injector');
W   = addWell(W, G.Matrix, G.Matrix.rock, cellprod, 'Type', 'bhp', ...
    'Val', 100*barsa, 'Sign', -1, 'Comp_i', [0, 1], 'Name', 'Producer');
plotWell(G,W);

%% Initialize state variables
% Once the transmissibilities are computed, we can generate the
% transmissiblity matrix 'A' using the 'two-point flux approximation'
% scheme and initialize the solution structure.
dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 0);
state.wellSol = initWellSol(W, 0);
[A,q] = getSystemIncompTPFA(state, G, T, fluid, 'use_trans', true); 

%% Setup multiscale grids
% Next, we define a 5-by-5-by-5 matrix coarse grid such that each coarse
% block contains 5-by-5-by-5 fine cells. The fracture plane is partitioned
% into 3-by-4 coarse blocks leading to a coarsening ratio of around 30.
% Additionally, we also define the support regions for the fracture and
% matrix basis functions. Fracture support region is defined based on a
% topological distance based algorithm. The matrix and fracture coarse
% grids are plotted in the next section.
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


%% Compute initial pressure
dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state_fs = incompTPFA(state, G, T, fluid,  ...
    'Wells', W, 'MatrixOutput', true, 'use_trans',true);

%% Compute basis functions
% Using the transmissibility matrix 'A' we compute the basis functions for
% each fracture and matrix coarse block using the restriction smoothed
% basis method. Note that the matrix 'A' does not contain any source terms
% or boundary conditions. They are added to the coarse linear system when
% computing the multiscale pressure in the next section.

dispif(mrstVerbose, 'Computing basis functions...\n\n');
basis_sb = getMultiscaleBasis(CG, A, 'type', 'msrsb');
clf; plotToolbar(G,basis_sb.B,'filterzero',true); view(-135,30)
plotGrid(CG,'FaceColor','none','LineWidth',1);
axis tight; colormap(jet); colorbar;
title('Basis Functions in the matrix');

%% Compute multiscale solution
dispif(mrstVerbose, 'Computing multiscale solution...\n\n');
[state_ms,~] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'Wells', W,'use_trans',true);

%% Plot initial pressure
figure;
plotToolbar(G, state_fs.pressure); plotWell(G,W); 
colormap jet; colorbar
view(15,20)
axis tight off
title('Fine scale')

figure;
plotToolbar(G, state_ms.pressure); plotWell(G,W); 
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
% We solve the two-phase system using a sequential splitting in which the
% pressure and fluxes are computed by solving the flow equation and then
% held fixed as the saturation is advanced according to the transport
% equation. Note that the flow equations are solved in the fine scale as
% well as on the coarse scale using an algebraic multiscale strategy. The
% multiscale flux field obtained at fine scale resolution is reconstructed
% to be conservative before solving the transport equation. This procedure
% is repeated for a given number of time steps (here we use 30 equally
% spaced time steps amounting to 90 % PV Injected). The error introduced by
% this splitting of flow and transport can be reduced by iterating each
% time step until e.g., the residual is below a certain user-prescribed
% threshold (this is not done herein).
pv     = poreVolume(G,G.rock);
nt     = 30;
Time   = 0.9*(sum(pv)/abs(sum(state_fs.wellSol(1).flux)));
dT     = Time/nt;

pvi = zeros(nt,1);
sol_fs = cell(nt,1); sol_ms = cell(nt,1);
e = zeros(nt,1);

% Prepare plotting
clf, set(gcf,'Position',[0 0 980 400]); colormap(flipud(winter));
txt = {'Reference','F-MsRSB'};
for i=1:2, subplot(1,2,i), 
   plotGrid(cartGrid([1 1 1],physdim), 'FaceColor','none','EdgeColor','k','LineWidth',1);
   plotWell(G,W); axis tight; view(20,15); title(txt{i});
end
[hfs,hms]=deal([]);

t  = 0;
B = basis_sb.B;
R = controlVolumeRestriction(CG.partition);
count = 1;
hwb = waitbar(0,'Time loop');
while t < Time,
    state_fs = implicitTransport(state_fs, G, dT, G.rock, fluid, 'wells', W, 'Trans', T,'verbose',true);
    state_ms = implicitTransport(state_ms, G, dT, G.rock, fluid, 'wells', W, 'Trans', T);
    % Check for inconsistent saturations
    s = [state_fs.s(:,1); state_ms.s(:,1)];
    assert(max(s) < 1+eps && min(s) > -eps);
    
    % Plot solution
    subplot(1,2,1); delete(hfs);
    hfs = plotCellData(G, state_fs.s(:,1),state_fs.s(:,1)>1e-3);
    subplot(1,2,2); delete(hms);
    hms = plotCellData(G, state_fs.s(:,1),state_ms.s(:,1)>1e-3);
   
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
    waitbar(t/Time, hwb);
    pvi(count) = 100*(sum(state_fs.wellSol(1).flux)*t)/sum(pv);
    e(count,1) = sum(abs(state_fs.s(:,1) - state_ms.s(:,1)).*pv)/sum(pv.*state_fs.s(:,1));
    
    fprintf([num2str(pvi(count)), '%% PV injected \n']);
    count = count + 1;
end
close(hwb);

%% Plot saturations
% To better see the evolution of the saturation front, you should change
% the patch properties (click the colored triangle in the bottom row of the
% menu) and set the faces to be semi-transparent (e.g., a value .5 to .7)
figure; plotToolbar(G,sol_fs); colormap(flipud(winter)); plotWell(G,W); view(15,20);
figure; plotToolbar(G,sol_ms); colormap(flipud(winter)); plotWell(G,W); view(15,20)

%% Plot error in saturation 

figure;
plot(pvi,e*100, '--+b');
ylabel('e [%]')
xlabel('PVI [%]'); 
set(gca,'XGrid','on','YGrid','on');
axis tight

e_eq = '$$ e = \frac{ \sum ( |S_w^{fs}-S_w^{f-msrsb}| \times pv) }{ \sum (S_w^{fs} \times pv) } $$';
title(e_eq,'interpreter','latex');

% <html>
% <p><font size="-1">
% Copyright 2009-2016 TU Delft and SINTEF ICT, Applied Mathematics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>