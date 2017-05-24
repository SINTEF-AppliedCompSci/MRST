%% F-MsRSB applied to a heterogeneous 2D domain
% Two-phase example modeling water injection through a quarter 5-spot into a
% 2-dimensional fractured porous media. The flow problem is solved both by a
% fine-scale and a multiscale solver.

% Notice that you need to have Metis installed to get this example to work.
% To get Metis working, you also need to set the global variable METISPATH.
% This can be done in your 'startup_user.m' file.

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add coarsegrid;      % functionality for coarse grids
mrstModule add ad-core;         % NNC support for coarse grids
mrstModule add msrsb;           % MsRSB solvers
mrstModule add mrst-gui;        % plotting routines
mrstModule add incomp;          % Incompressible fluid models
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path

%% Grid and fracture line(s)
% Construct a Cartesian grid comprising 100-by-100 cells, where each cell
% has dimension 1-by-1 m^2. Define 1 fracture line, roughly 85 m in length,
% extending diagonally through the centre of the domain.

celldim = [100 100];
physdim = [100 100];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

fl = [20 20 80 80];

%% Process fracture line(s)
% Using the input fracture lines, we identify independent fracture networks
% comprising of connected lines. In this example, there is only 1 fracture
% network consisting of 1 fracture line. We also identify the fine-cells
% in the matrix containing these fractures. Fracture aperture is set to
% 0.04 meters. The matrix grid and fracture lines are plotted.

dispif(mrstVerbose, 'Processing user input...\n\n');
[G,fracture] = processFracture2D(G,fl);
fracture.aperture = 1/25; % Fracture aperture

figure;
plotFractureLines(G,fracture);
box on

%% Compute CI and construct fracture grid
% For each matrix block containing a fracture, we compute a fracture-matrix
% conductivity index (CI) in accordance with the hierarchical fracture
% model (a.k.a. embedded discrete fracture model). Following that, we
% compute a fracture grid where the fracture cell size is defined to be
% 1 m. Next, the fracture grid is plotted on top of the matrix grid.

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 0.9; cell_size = 1; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); box on

%% Set rock properties in fracture and matrix
% For this example, we will extract permeability and porosity the first
% layer of the SPE10 dataset and re-sample the same for the matrix grid at
% the appropriate grid resolution. Fracture permeability is set to 1000
% Darcy with 50% porosity in each fracture grid cell. The rock properties
% are then plotted.

dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');

mrstModule add spe10
rock = SPE10_rock(1);
p = reshape(rock.perm(:,1)*milli*darcy,60,[]);
perm = sampleFromBox(G,p');
p = reshape(rock.poro,60,[]);
poro = sampleFromBox(G,p');
% Since we don't use actnum to denote N/G or active cells we ignore cells
% with very low porosity.
poro(poro<=0.1) = 0.1; 
G.rock = makeRock(G,perm(:),poro(:));
K_frac = 1000*darcy; % Darcy
poro_frac = 0.5;
for i = 1:numel(fieldnames(G.FracGrid))
    G.FracGrid.(['Frac',num2str(i)]).rock.perm = K_frac*ones(G.FracGrid.(['Frac',num2str(i)]).cells.num,1);
    G.FracGrid.(['Frac',num2str(i)]).rock.poro = poro_frac*ones(G.FracGrid.(['Frac',num2str(i)]).cells.num,1);
end
clf; plotToolbar(G,G.rock); colormap(jet); colorbar

%% Define fluid properties
% Define a two-phase fluid model without capillarity. The fluid model has
% the values:
%
% * densities: [rho_w, rho_o] = [1000 700] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 2] cP.
% * corey-coefficient: [2, 2] = [2 2].

fluid = initSimpleFluid('mu' , [   1,  2] .* centi*poise     , ...
    'rho', [1000, 700] .* kilogram/meter^3, ...
    'n'  , [   2,   2]);

%% Define fracture connections as NNC and compute the transmissibilities
% In this section, we use the function defineNNCandTrans to combine the
% fracture and matrix grid structures into a single grid structure. In
% addition to that, we assign a 'non-neighbouring connection (NNC)' status
% to every fracture-matrix connection. For example, if fracture cell 'i' is
% located in the matrix cell 'j', then the connection i-j is defined as an
% NNC. To compute the flux between these elements, we compute a
% transmissibility for each NNC using the CI's computed earlier. Vector T
% contains the transmissibility for each face in the combined grid and each
% NNC.

[G,T] = defineNNCandTrans(G,F,fracture);

%% Add wells
% We will include two wells, both the injector and the producer controlled
% by bottom-hole pressure. The injector and producer are located at the
% bottom left and top right corner of the grid, respectively. Wells are
% described using a Peaceman model, giving an extra set of equations that
% need to be assembled, see simpleWellExample.m for more details on the
% Peaceman well model.

W = addWell([], G.Matrix, G.Matrix.rock, 1,'InnerProduct', 'ip_tpf','Type', ...
    'bhp', 'Val', 120*barsa,'Radius', .1, 'Comp_i', [1, 0]);
W = addWell(W, G.Matrix, G.Matrix.rock, prod(celldim), 'InnerProduct', 'ip_tpf', 'Type', ...
    'bhp' , 'Val', 100*barsa, 'Radius', .1, 'Comp_i', [0, 1]);

%% Initialize state variables
% Once the transmissibilities are computed, we can generate the
% transmissiblity matrix 'A' using the 'two-point flux approximation'
% scheme and initialize the solution structure.

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 0);
state.wellSol = initWellSol(W, 0);
[A,q] = getSystemIncompTPFA(state, G, T, fluid, 'use_trans', true); 

%% Setup multiscale grids 
% Next, we define a 10-by-10 matrix coarse grid and 2 coarse blocks (or
% coarse degrees of freedom) in the fracture. Additionally, we also define
% the support regions for the fracture and matrix basis functions. The
% matrix and fracture coarse grids are plotted.

dispif(mrstVerbose, 'Defining coarse grids and interaction regions...\n\n');
coarseDims = [10 10]; % coarsening factor in each direction
dof_frac = 2; % Fracture dof per fracture network
[CG, CGf] = getRsbGridsHFM(G, fracture.network, 'coarseDims', coarseDims,...
    'dof_frac',dof_frac, 'paddedPartition', true, ...
    'coarseNodeOption', 'useCoarseCellEndPoints'); 
clf; plotFractureCoarseGrid2D(G,CG.partition,F)

%% Compute initial pressure

dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state_fs = incompTPFA(state, G, T, fluid,  ...
    'Wells',W, 'MatrixOutput', true, 'use_trans',true);

%% Compute basis functions
% Using the transmissibility matrix 'A' we compute the basis functions for
% each fracture and matrix coarse block using the restriction smoothed
% basis method. Note that the matrix 'A' does not contain any source terms
% or boundary conditions. They are added to the coarse linear system when
% computing the multiscale pressure in the next section.

dispif(mrstVerbose, 'Computing basis functions...\n\n');
basis_sb = getMultiscaleBasis(CG, A, 'type', 'rsb');
clf; plotToolbar(G,basis_sb.B,'filterzero',true);
prm = log10(G.rock.perm); mx = max(prm); mn = min(prm);
G.nodes.z = 1e-3*ones(G.nodes.num,1);
plotCellData(G,(prm-mn)./(mx-mn),'EdgeColor','none','FaceAlpha',.4);
G.nodes = rmfield(G.nodes,'z');
plotGrid(CG,'FaceColor','none');
axis tight; colorbar;
title('Basis functions plotted in the matrix');

%% Compute multiscale solution

dispif(mrstVerbose, 'Computing multiscale solution...\n\n');
[state_ms,~] = incompMultiscale(state, CG, T, fluid, basis_sb,...
    'Wells', W,'use_trans',true);

%% Incompressible Two-Phase Flow
% We solve the two-phase system using a sequential splitting in which the
% pressure and fluxes are computed by solving the flow equation and then
% held fixed as the saturation is advanced according to the transport
% equation. Note that the flow equations are solved in the fine scale as
% well as on the coarse scale using an algebraic multiscale strategy. The
% multiscale flux field obtained at fine scale resolution is reconstructed
% to be conservative before solving the transport equation. This procedure
% is repeated for a given number of time steps (here we use 45 equally
% spaced time steps). The error introduced by this splitting of flow and
% transport can be reduced by iterating each time step until e.g., the
% residual is below a certain user-prescribed threshold (this is not done
% herein).

pv     = poreVolume(G,G.rock);
nt     = 45;
Time   = 0.5*(sum(pv)/state_fs.wellSol(1).flux);
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
figure; set(gcf,'Position',[0 0 800 400]); colormap(flipud(winter));
for i=1:2, subplot(1,2,i),
   plotCellData(G,(prm-mn)./(mx-mn),'EdgeColor','none','FaceAlpha',.4);
   drawnow;
end
[hms,hfs]=deal([]);
hwb = waitbar(0,'Time loop');
while t < Time,
    state_fs = implicitTransport(state_fs, G, dT, G.rock, fluid, 'wells', W, 'Trans', T,'verbose',true);
    state_ms = implicitTransport(state_ms, G, dT, G.rock, fluid, 'wells', W, 'Trans', T);
    % Check for inconsistent saturations
    s = [state_fs.s(:,1); state_ms.s(:,1)];
    assert(max(s) < 1+eps && min(s) > -eps);

    % Plot solutions
    subplot(1,2,1);
    delete(hfs);
    hfs=plotCellData(G,state_fs.s(:,1),state_fs.s(:,1)>0,'EdgeColor','none');
    line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
    view(0,90); caxis([0 1]); title('Fine scale');
    
    subplot(1,2,2);
    delete(hms);
    hms=plotCellData(G,state_ms.s(:,1),state_ms.s(:,1)>0,'EdgeColor','none');
    line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
    view(0,90); caxis([0 1]); title('F-MsRSB');
    drawnow

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
    waitbar(t/Time,hwb);
    pvi(count) = 100*(sum(state_fs.wellSol(1).flux)*t)/sum(pv);
    e(count,1) = sum(abs(state_fs.s(:,1) - state_ms.s(:,1)).*pv)/sum(pv.*state_fs.s(:,1));
    pfs(count,1) = state_fs.s(W(2).cells,1);
    pms(count,1) = state_ms.s(W(2).cells,1);
    count = count + 1;
    
end
close(hwb);

%% Plot saturations

plotNo = 1;
figure; hold on; colormap(flipud(winter))
boundary = any(G.Matrix.faces.neighbors==0,2);
facelist = 1:G.Matrix.faces.num;
bfaces = facelist(boundary);
for i = nt/3:nt/3:nt
    state_fs = sol_fs{i,1}; state_ms = sol_ms{i,1};
    heading = [num2str(round(pvi(i))),  ' % PVI'];
    % Plot saturation
    r = 0.01;
    subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.48]), cla
    plotCellData(G,state_fs.s(:,1),'EdgeColor','none');
    line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
    plotFaces(G,bfaces,'k','linewidth',1)
    axis square off, 
    title(['Reference: ', heading],'FontSize',8);
    view(0,90); caxis([0 1]);
    
    subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.48]), cla
    plotCellData(G,state_ms.s(:,1),'EdgeColor','none');
    line(fl(:,1:2:3)',fl(:,2:2:4)','Color','r','LineWidth',0.5);
    plotFaces(G,bfaces,'k','linewidth',1)
    axis square off, 
    title(['F-MsRSB: ',  heading],'FontSize',8);
    view(0,90); caxis([0 1]);
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