mrstModule add dg vem vemmech ad-props ad-core ad-blackoil       ...
    blackoil-sequential gasinjection mrst-gui reorder matlab_bgl ...
    coarsegrid incomp upr msrsb
mrstVerbose on;

%%

n  = 50;
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
p_cart = partitionUI(G_cart, [20, 20]);
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
modelDG  = getSequentialModelFromFI(modelFI);
modelASI = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, modelSIF.transportModel, G);
% modelASI.transportModel = TransportOilWaterModel(G, rock, fluid);


if 0
    s = getSmootherFunction('type', 'ilu', 'iterations', 1);

    mssolver = MultiscaleVolumeSolverAD(G, 'getSmoother', s);
    msmodel = MultiscalePressureModel(modelASI.pressureModel.G,...
                                      modelASI.pressureModel.rock,...
                                      modelASI.pressureModel.fluid,...
                                      modelASI.pressureModel, mssolver);
    modelASI.pressureModel = msmodel;
end

[jt, ot, mt] = deal(Inf);
mt = 1e-3;
ot = 1e-3;
degree = 1;
if 0
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
transportModel = TransportOilWaterModelDG(GF, rockF, fluid, ...
                                   'disc'    , disc        , ...
                                   'dsMaxAbs', 0.2, ...
                                   'nonlinearTolerance', 1e-3);

modelASI = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, transportModel, G);
end                          
                               
discDG = DGDiscretization(modelDG.transportModel               , ...
                        'degree'               , degree       , ...
                        'basis'                , 'legendre'   , ...
                        'useUnstructCubature'  , true        , ...
                        'jumpTolerance'        , jt           , ...
                        'outTolerance'         , ot           , ...
                        'outLimiter'           , 'kill'       , ...
                        'meanTolerance'        , mt           , ...
                        'limitAfterConvergence', false        , ...
                        'plotLimiterProgress'  , false        );
modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, ...
                                   'disc'    , discDG        , ...
                                   'dsMaxAbs', 0.2, ...
                                   'nonlinearTolerance', 1e-3);
                               


% pModel = modelASI.pressureModel;
% msSolver = MultiscaleVolumeSolverAD(G);
% modelASI.pressureModel = MultiscalePressureModel(GF, rockF, fluid, pModel, msSolver);
% modelASI.pressureNonLinearSolver.LinearSolver = msSolver;

%%

time = 2*year;
rate = 2*sum(poreVolume(GF, rockF))/time;

xw = [0,0; l,l];
% xw = G.cells.centroids([1, G.cells.num], :);

W = [];
d = pdist2(G.cells.centroids, xw);

[~, ix] = min(d(:,1));
W = addWell(W, G, rock, ix, 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
[~, ix] = min(d(:,2));
W = addWell(W, G, rock, ix, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

dt    = 7*day;
dtvec = rampupTimesteps(time, dt, 0);
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

sW      = 0.0;
state0F = initResSol(GF, 100*barsa, [sW,1-sW]);
state0F.bfactor = [fluid.bW(state0F.pressure), fluid.bO(state0F.pressure)];

%%

% modelASI.storeGrids = true;
% modelASI.plotProgress = true;
% modelASI.pressureModel.extraStateOutput = true;
% 
% state0F.transportState        = state0;
% state0F.transportState.degree = degree*ones(G.cells.num, 1);
% state0F.transportState.G      = G;
% state0F.transportState.pv     = modelASI.transportModel.operators.pv;
% state0F.transportState        = disc.assignDofFromState(state0F.transportState);
% state0F.transportState        = disc.updateDofPos(state0F.transportState);
% state0F.transportModel        = modelASI.transportModel;
% % 
% state0 = disc.assignDofFromState(state0);
% state0 = disc.updateDofPos(state0);
% state0F.sdof = state0.sdof;
% % state0F.nDof = state0.nDof;


modelASI.storeGrids                     = true;
modelASI.plotProgress                   = true;
modelASI.pressureModel.extraStateOutput = true;

state0F.transportState        = state0;
% state0F.transportState.degree = degree*ones(G.cells.num, 1);
state0F.transportState.G      = G;
state0F.transportState.pv     = modelASI.transportModel.operators.pv;
% state0F.transportState        = disc.assignDofFromState(state0F.transportState);
% state0F.transportState        = disc.updateDofPos(state0F.transportState);
state0F.transportModel        = modelASI.transportModel;


[wsASI, stASI, rep] = simulateScheduleAD(state0F, modelASI, scheduleF);

%%

[wsSI, stSI, rep] = simulateScheduleAD(state0, modelSI, schedule);

%%

state0 = assignDofFromState(modelDG.transportModel.disc, state0);
[wsDG, stDG, rep] = simulateScheduleAD(state0, modelDG, schedule);

%%

[wsSIF, stSIF, rep] = simulateScheduleAD(state0F, modelSIF, scheduleF);

%%

ws = {wsSI, wsDG, wsASI, wsSIF};
names    = {'Coarse', ['Coarse dG(', num2str(degree), ')'], ...
            ['Adaptive dG(', num2str(degree), ')'], 'Fine'};
        
pIx = [2];
ws = ws(pIx);
names = names(pIx);
        
plotWellSols(ws, 'datasetNames', names);

%%

close all
plotToolbar(GF, stASI)
axis equal tight
colormap(jet);

%%

close all
plotToolbar(G, stDG)
axis equal tight

%%

close all
plotToolbar(G, stSI)
axis equal tight
colormap(jet);

%%

wcut = nan(numel(dt,1),3);
for sNo = 1:numel(stASI)
    wcut(sNo,:) = [wsSI{sNo}(2).wcut, wsASI{sNo}(2).wcut, wsSIF{sNo}(2).wcut];
end

%%



close all
fig = figure('Position', [-2000, 0, 1500, 600]);
M = struct('cdata',[],'colormap',[]);
d = 10;

clr = lines(3);
    
dt = cumsum(schedule.step.val)/day;

% subplot(1,3,3)
% hold on
% wc = wcut;
% wc(2:end,:) = nan;
% xlim([0, dt(end)]);
% ylim([0,1]);
% hold off
% pbaspect([1,1,1])
% box on
% ylabel('Water cut')
% xlabel('Time (days)')
% legend({'Coarse', 'Adaptive', 'Reference'}, 'location', 'northwest');

for sNo = 1:numel(stASI)
    
    subplot(1,3,1)
    cla;
    hold on
    plotCellData(GF, stASI{sNo}.s(:,1), 'edgec', 'none')
    plotGrid(stASI{sNo}.G, 'facecolor', 'none')
    hold off
    caxis([0,1]);
    axis equal tight
    
    subplot(1,3,2)
    plotCellData(G, stSI{sNo}.s(:,1), 'edgec', 'none')
    plotGrid(G, 'facec', 'none', 'edgealpha', 0.2);
    caxis([0,1]);
    axis equal tight
    
    colormap(jet)
    
    subplot(1,3,3)
    wc = wcut;
    wc(sNo+1:end,:) = nan;
    hold on
    plot(dt, wc(:,1), 'linew', 2, 'color', clr(1,:));
    plot(dt, wc(:,2), 'linew', 2, 'color', clr(2,:));
    plot(dt, wc(:,3), '--', 'linew', 4, 'color', clr(3,:));
    hold off
    pbaspect([1,1,1])
    xlim([0, dt(end)]);
    ylim([0,1]);
    box on
    ylabel('Water cut')
    xlabel('Time (days)')
    legend({'Coarse', 'Adaptive', 'Reference'}, 'location', 'northwest');
%     for wNo = 1:3
%         h(wNo).YData = wc(:,wNo);
%     end

    if 0
        rect = [d, d, fig.Position(3:4) - [d,d]];
        M(sNo) = getframe(fig, rect);
    end
    
    pause(0.05)
    
end

%%



pth = mrstPath('dg');
name = 'refinement';
duration = 10;
vo = VideoWriter(fullfile(pth, name));
vo.FrameRate = numel(stASI)/duration;
open(vo);

writeVideo(vo, M);

close(vo)

%%

plotWellSols({wsSI, wsASI, wsSIF})

%%

close all
plotToolbar(G, st)