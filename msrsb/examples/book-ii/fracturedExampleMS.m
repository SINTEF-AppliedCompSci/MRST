%% Injection of a highly mobile fluid into a fractured medium
% This example is a simple incompressible, two-phase analogue of the
% examples studied in Section 5.2 of Møyner & Tchelepi, SPE J, 23(6), 2018,
% doi: 10.2118/182679-PA. We inject a fluid of density 300 kg/m^3 and
% viscosity 0.3 cP into another fluid of density 1000 kg/m^3 and viscosity
% 1 cP that fills a layered medium containing thirteen high-permeability
% fractures.
mrstModule add coarsegrid incomp msrsb book ad-core

%%  Set up the fine-scale model
gravity reset off
rdim = [1000 500];
load setup_fracture.mat;
G.nodes.coords = bsxfun(@times, G.nodes.coords, rdim);
G = computeGeometry(G);
perm(G.cells.tag)=5000;
rock = makeRock(G, perm*milli*darcy, 0.3);
hT   = computeTrans(G, rock);

figure
plotCellData(G, log10(perm),'EdgeAlpha',.1);
mrstColorbar(perm,'south',true); axis tight

fluid = initSimpleFluid('mu', [.3, 1]*centi*poise, ...
    'n', [2, 2], 'rho', [300, 1000]*kilo*gram/meter.^3);

state0 = initResSol(G, 300*barsa, [0, 1]);

icell = findEnclosingCell(G, [0.1, 0.1].*rdim);
pcell = findEnclosingCell(G, [0.9, 0.9].*rdim);
pv    = sum(poreVolume(G, rock));
time  = 10*year;

W  = addWell([], G, rock, icell, 'type', 'rate', ...
             'val',  pv/time, 'compi', [1 0], 'Name', 'I');
W  = addWell(W, G, rock, pcell, 'type', 'rate', ...
             'val', -pv/time, 'compi', [1 0], 'Name', 'P');
plotGrid(G,[icell pcell],'FaceColor','w');

%% Set up multiscale problem
ps    = sampleFromBox(G, reshape(1:8*7,[7 8]));
p     = processPartition(G, compressPartition(ps));
%pf    = zeros(G.cells.num,1); pf(G.cells.tag>0)=1;
%p     = processPartition(G, compressPartition(pf*8*7 + ps));
CG    = generateCoarseGrid(G, compressPartition(p));
CG    = coarsenGeometry(CG);
CG    = storeInteractionRegion(CG);
basis = getMultiscaleBasis(CG, getIncomp1PhMatrix(G, hT), 'type', 'msrsb');

%% Set up solvers
dt = rampupTimesteps(time, time/100, 5);

psolver  = @(u)     incompTPFA(u, G, hT, fluid, 'W', W);
tsolver  = @(u, dt) implicitTransport(u, G, dt, rock, fluid, 'W', W);
mssolve  = @(u)     incompMultiscale(u, CG, hT, fluid, basis, 'W', W);

[rstates, rws, mstates, mws] = deal(cell(numel(dt)+1,1));
rstates{1} = psolver(state0); rws{1} = getWellSol(W, rstates{1}, fluid);
mstates{1} = mssolve(state0); mws{1} = getWellSol(W, mstates{1}, fluid);

%% Plot the initial saturation distribution
figure('position', [580, 540, 950, 230])
subplot(1,2,1)
hp(1) = plotCellData(G, rstates{1}.s(:,1),'EdgeAlpha',.1); 
plotGrid(G,[icell pcell],'FaceColor','w');
axis tight, caxis([0 1]), colormap(flipud(winter))
title('Fine scale','FontSize',12,'FontWeight','normal');

subplot(1,2,2)
hp(2) = plotCellData(G, mstates{1}.s(:,1),'EdgeAlpha',.1); 
axis tight, caxis([0 1]), colormap(flipud(winter(10).^1.5))
plotGrid(G,[icell pcell],'FaceColor','w');
plotFaces(CG,1:CG.faces.num,'EdgeColor','k','FaceColor','none');
title('MsRSB: no iterations','FontSize',12,'FontWeight','normal');

%% Run the simulation
for i = 1:numel(dt)
    state        = psolver(rstates{i});
    rstates{i+1} = tsolver(state, dt(i));
    rws{i+1}     = getWellSol(W, rstates{i+1}, fluid);
    hp(1).CData  = rstates{i+1}.s(:,1);
    
    state        = mssolve(mstates{i});
    mstates{i+1} = tsolver(state, dt(i));
    mws{i+1}     = getWellSol(W, mstates{i+1}, fluid);
    hp(2).CData  = mstates{i+1}.s(:,1);
    
    drawnow
end
plotWellSols({rws, mws}, [0; cumsum(dt)], 'datasetnames', {'fine scale','msrsb'});