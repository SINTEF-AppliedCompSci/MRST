%% Simulate gravity segregation with the MsRSB method
% This example is a simple demonstration in how to use the multiscale
% solver to solve a gravity segregation problem using the MsRSB method.

mrstModule add coarsegrid incomp msrsb

%%  Set up the fine-scale model
rng(07521362);
gravity reset on
G    = computeGeometry(cartGrid([20 1 40],[500 20 60]*meter));
p    = gaussianField(G.cartDims, [0.2 0.4], [11 1 3], 2.5);
K    = p.^3.*(1e-5)^2./(0.81*72*(1-p).^2);
rock = makeRock(G, K(:), p(:));
hT   = computeTrans(G, rock);

figure
K = convertTo(rock.perm(:,1),milli*darcy);
plotCellData(G, K,'EdgeAlpha',.1); mrstColorbar(K,'south');
axis tight, view(0,0)

fluid = initSimpleFluid('mu', [1, .5]*centi*poise, ...
    'n', [2, 2], 'rho', [1000, 600]*kilo*gram/meter.^3);

state0 = initResSol(G, 0, [0, 1]);
ind    = G.cells.centroids(:,3)<30;
state0.s(ind,:) = 1 - state0.s(ind,:);


%% Set up multiscale problem
coarsen = [5 1 5];
coarsedims = ceil(G.cartDims./coarsen);
p     = partitionUI(G, ceil(G.cartDims./coarsen));
CG    = generateCoarseGrid(G, compressPartition(p));
CG    = coarsenGeometry(CG);
CG    = storeInteractionRegionCart(CG);
basis = getMultiscaleBasis(CG, getIncomp1PhMatrix(G, hT), 'type', 'msrsb');
smoother = getSmootherFunction('type', 'ilu0','iterations',1);
iargs = {'tolerance', 1e-3, 'iterations', 50,...
          'useGMRES', true, 'getSmoother', smoother};

%% Set up solvers
time = 10*year;
Nt   = 100;
dt   = time/Nt;

psolver  = @(state) incompTPFA(state, G, hT, fluid);
tsolver  = @(state) implicitTransport(state, G, dt, rock, fluid);
mssolve  = @(state) incompMultiscale(state, CG, hT, fluid, basis);
isolver  = @(state) incompMultiscale(state, CG, hT, fluid, basis, iargs{:});

[rstates,mstates,istates] = deal(cell(Nt+1,1));
rstates{1} = psolver(state0);
mstates{1} = mssolve(state0);
istates{1} = isolver(state0);

%% Plot the initial saturation distribution
figure('position', [50, 50, 1000, 500])
subplot(1,3,1)
hp(1) = plotCellData(G, rstates{1}.s(:,1),'EdgeAlpha',.1); 
view(0,0); axis tight, caxis([0 1]), colormap(flipud(winter))
title('Fine scale');

subplot(1,3,2)
hp(2) = plotCellData(G, mstates{1}.s(:,1),'EdgeAlpha',.1); 
view(0,0); axis tight, caxis([0 1]), colormap(flipud(winter(10).^1.5))
plotFaces(CG,1:CG.faces.num,'EdgeColor','k','FaceColor','none');
title('MsRSB: no iterations');

subplot(1,3,3)
hp(3) = plotCellData(G, istates{1}.s(:,1),'EdgeAlpha',.1); 
view(0,0); axis tight, caxis([0 1]), colormap(flipud(winter(10).^1.5))
plotFaces(CG,1:CG.faces.num,'EdgeColor','k','FaceColor','none');
title('MsRSB: iterative GMRES');
[~,c] = boundaryFaces(G);

%% Run the simulation
for i = 1:Nt    
    state        = psolver(rstates{i});
    rstates{i+1} = tsolver(state);
    hp(1).CData  = rstates{i+1}.s(c,1);
    
    state        = mssolve(mstates{i});
    mstates{i+1} = tsolver(state);
    hp(2).CData  = mstates{i+1}.s(c,1);
    
    state        = isolver(istates{i});
    istates{i+1} = tsolver(state);
    hp(3).CData  = istates{i+1}.s(c,1);
    drawnow
end