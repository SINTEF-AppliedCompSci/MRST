mrstModule add dg vem vemmech ad-props ad-core ad-blackoil ...
    blackoil-sequential gasinjection reorder matlab_bgl upr coarsegrid ...
    mrst-gui mrst-utils

mrstVerbose on

%%

n  = 30;
l = 1000;
GF = pebiGrid(l/n, [l,l]);
close all
plotGrid(GF)
axis equal tight

%%

GF = computeGeometry(GF);
GF = computeCellDimensions2(GF);

% Cartesian coarse grid
G_cart = cartGrid([100, 100]);
p_cart = partitionUI(G_cart, [10, 10]);
p_cart = sampleFromBox(GF, reshape(p_cart, G_cart.cartDims));
GC = generateCoarseGrid(GF, p_cart);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC);
% GC = storeInteractionRegionCart(GC);

%%

close all

mappings = getRefinementMappings(GC, GC, GF, [1, GC.cells.num]);
G        = generateCoarseGrid(GF, mappings.newPartition);
G        = coarsenGeometry(G);
G.cells.refined = mappings.refined;

% [G, map1] = refineGrid(GC, GC, GF, [1, GC.cells.num]');
G = coarsenCellDimensions(G);

G.cells.ghost = false(G.cells.num,1);
GF.cells.ghost = false(GF.cells.num,1);


%%

gravity reset off

rock  = makeRock(G, 100*milli*darcy, 1);
rockF = makeRock(GF, 100*milli*darcy, 1);

fluid = initSimpleADIFluid('phases', 'WO'                      , ...
                           'rho'   , [1000, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise    , ...
                           'n'     , [1, 1]                    );

modelFIF = TwoPhaseOilWaterModel(GF, rockF, fluid);
modelSIF = getSequentialModelFromFI(modelFIF);
modelFI  = TwoPhaseOilWaterModel(G, rock, fluid);
modelSI  = getSequentialModelFromFI(modelFI);
modelASI = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, modelSIF.transportModel, G);

transportModel = modelSIF.transportModel;

[jt, ot, mt] = deal(Inf);
mt = 1e-3;
ot = 1e-3;
degree = 1;
disc = DGDiscretization(modelASI.transportModel               , ...
                        'degree'               , degree       , ...
                        'basis'                , 'legendre'   , ...
                        'useUnstructCubature'  , false        , ...
                        'jumpTolerance'        , jt           , ...
                        'outTolerance'         , ot           , ...
                        'outLimiter'           , 'orderReduce', ...
                        'meanTolerance'        , mt           , ...
                        'limitAfterConvergence', false        , ...
                        'plotLimiterProgress'  , false        );
transportModelDG = TransportOilWaterModelDG(GF, rockF, fluid, ...
                                   'disc'    , disc        , ...
                                   'dsMaxAbs', 0.2, ...
                                   'nonlinearTolerance', 1e-3);

modelASIDG = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, transportModelDG, G);

tm = ReorderingModelDG(transportModelDG);
modelASIDGreorder = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, tm, G);

% modelASIreorder = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, transportModel, G);
modelASIDGreorder.transportModel.chunkSize = 100;
modelASIDGreorder.transportModel.parent.extraStateOutput = true;

%%

time = 1*year;
rate = 1*sum(poreVolume(GF, rockF))/time;

xw = [0,0; l,l];
% xw = G.cells.centroids([1, G.cells.num], :);

W = [];
d = pdist2(G.cells.centroids, xw);

[~, ix] = min(d(:,1));
W = addWell(W, G, rock, ix, 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
[~, ix] = min(d(:,2));
W = addWell(W, G, rock, ix, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

dt    = 7*day;
dtvec = rampupTimesteps2(time, dt, 0);
% dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);

WF = [];
d = pdist2(GF.cells.centroids, xw);

[~, ix] = min(d(:,1));
WF = addWell(WF, GF, rockF, ix, 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
[~, ix] = min(d(:,2));
WF = addWell(WF, GF, rockF, ix, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

scheduleF = simpleSchedule(dtvec, 'W', WF);

state0F = initResSol(GF, 100*barsa, [sW,1-sW]);
state0F.bfactor = [fluid.bW(state0F.pressure), fluid.bO(state0F.pressure)];

%%

modelASI.storeGrids                     = true;
modelASI.plotProgress                   = true;
modelASI.pressureModel.extraStateOutput = true;

state0F.transportState        = state0;
state0F.transportState.degree = degree*ones(G.cells.num, 1);
state0F.transportState.G      = G;
state0F.transportState.pv     = modelASI.transportModel.operators.pv;
state0F.transportState        = disc.assignDofFromState(state0F.transportState);
state0F.transportState        = disc.updateDofPos(state0F.transportState);
state0F.transportModel        = modelASI.transportModel;

%%

[wsASI, stASI, rep] = simulateScheduleAD(state0F, modelASI, scheduleF);

%%

state0F.transportModel  = modelASIDG.transportModel;
state0F.transportState.G = G;
[wsASIDG, stASIDG, rep] = simulateScheduleAD(state0F, modelASIDG, scheduleF);

%%

modelASIDGreorder.transportModel.parent.extraStateOutput = true;
modelASIDGreorder.pressureModel.extraStateOutput = true;
modelASIDGreorder.transportModel.plotProgress = true;
modelASIDGreorder.transportModel.plotAfterTimestep = false;
modelASIDGreorder.storeGrids = true;
modelASIDGreorder.plotProgress = false;
state0F.transportModel        = modelASIDGreorder.transportModel;
[wsASIreorder, stASIreorder, rep] = simulateScheduleAD(state0F, modelASIDGreorder, scheduleF);

%%

[wsSI, stSI, rep] = simulateScheduleAD(state0, modelSI, schedule);

%%

[wsSIF, stSIF, rep] = simulateScheduleAD(state0F, modelSIF, scheduleF);

%%

close all
cmap = mrstColormap('type', 'wateroil');

%%

ws = {wsSI, wsASIreorder, wsASI, wsASIDG, wsSIF};
names    = {'Coarse', ...
            ['Adaptive dG(', num2str(degree), ') reorder'], ...
            ['Adaptive dG(0)']       , ...
            ['Adaptive dG(', num2str(degree), ')']       , ...
            'Fine'};
        
pIx = 1:4;
ws = ws(pIx);
names = names(pIx);
        
plotWellSols(ws, scheduleF.step.val, 'datasetNames', names);

%%

figure
plotToolbar(GF, stASI)
axis equal tight
colormap(cmap);

%%

figure
plotToolbar(GF, stASIreorder)
axis equal tight
colormap(cmap);

%%



%%

close all
plotToolbar(G, stDG)
axis equal tight
colormap(jet);

%%

close all
plotToolbar(G, stSI)
axis equal tight
colormap(jet);

%%