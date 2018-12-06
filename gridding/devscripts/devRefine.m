mrstModule add dg vem vemmech ad-props ad-core ad-blackoil       ...
    blackoil-sequential gasinjection mrst-gui reorder matlab_bgl ...
    coarsegrid incomp
mrstVerbose on;

%%

n  = 100;
l = 1000;
GF = computeGeometry(cartGrid([n,n], [l,l]));
GF = computeCellDimensions2(GF);
rockF = makeRock(GF, 100*milli*darcy, 1);
p = partitionUI(GF, [10,10]);

GC = generateCoarseGrid(GF, p);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC);
GC = storeInteractionRegionCart(GC);

%%

% close all
[G, map1] = refineGrid(GC, GC, GF, [1, GC.cells.num]');
G.cells.ghost = false(G.cells.num,1);
GF.cells.ghost = false(GF.cells.num,1);
% plotGrid(G)
% G = GC;


%%

rock  = makeRock(G, 100*milli*darcy, 1);


fluid = initSimpleADIFluid('phases', 'WO'                      , ...
                           'rho'   , [1000, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise    , ...
                           'n'     , [2, 1]                    );

modelFIF  = TwoPhaseOilWaterModel(GF, rockF, fluid);
modelSIF  = getSequentialModelFromFI(modelFIF);                       
modelFI  = TwoPhaseOilWaterModel(G, rock, fluid);
modelSI  = getSequentialModelFromFI(modelFI);
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
modelASI.transportModel = TransportOilWaterModelDG(GF, rockF, fluid, ...
                                   'disc'    , disc        , ...
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

dt    = 30*day;
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

modelASI.pressureModel.extraStateOutput = true;
[wsASI, stASI, rep] = simulateScheduleAD(state0F, modelASI, scheduleF);

%%

[wsSI, stSI, rep] = simulateScheduleAD(state0, modelSI, schedule);

%%

[wsSIF, stSIF, rep] = simulateScheduleAD(state0F, modelSIF, scheduleF);

%%

plotWellSols({wsSI, wsASI, wsSIF});

%%

close all
plotToolbar(GF, stASI)

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
    hold on
    wc = wcut;
    wc(sNo+1:end,:) = nan;
    plot(dt, wc(:,1), 'linew', 2, 'color', clr(1,:));
    plot(dt, wc(:,2), 'linew', 2, 'color', clr(2,:));
    plot(dt, wc(:,3), '--', 'linew', 4, 'color', clr(3,:));
    xlim([0, dt(end)]);
    ylim([0,1]);
    hold off
    pbaspect([1,1,1])
    box on
    ylabel('Water cut')
    xlabel('Time (days)')
    legend({'Coarse', 'Adaptive', 'Reference'}, 'location', 'northwest');

    if 0
        rect = [d, d, fig.Position(3:4) - [d,d]];
        M(sNo) = getframe(fig, rect);
    end
    
    pause(0.1)
    
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