mrstModule add dg vem vemmech ad-props ad-core ad-blackoil       ...
    blackoil-sequential gasinjection mrst-gui reorder matlab_bgl ...
    coarsegrid incomp
mrstVerbose on;

%%

n  = 100;
l = 1000;
GF = computeGeometry(cartGrid([n,n], [l,l]));
% GF = computeGeometry(pebiGrid(l/n, [l,l]));
rockF = makeRock(GF, 100*milli*darcy, 1);

% T = getFaceTransmissibility(GF, rockF);
% p = partitionMETIS(GF, T, 25);
% p = processPartition(GF, p);
% p = compressPartition(p);
p = partitionUI(GF, [10,10]);

GC = generateCoarseGrid(GF, p);
GC = coarsenGeometry(GC);
GC = assignCoarseNodes(GC);

%%

[G, map1] = refineGrid(GC, GC, GF, [1, GC.cells.num]');
plotGrid(G)

% close all
% [G1, map1] = refineGrid(GC, GF, 1');
% figure, plotGrid(G1);
% 
% [G2, map2] = refineGrid(G1, GF, 3');
% figure, plotGrid(G2);
% G = G2;

% axis equal tight
% G = computeVEMGeometry(G)
%%

rock  = makeRock(G, 100*milli*darcy, 1);


fluid = initSimpleADIFluid('phases', 'WO'                      , ...
                           'rho'   , [1000, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise    , ...
                           'n'     , [1, 1]                    );

modelFI  = TwoPhaseOilWaterModel(G, rock, fluid);
modelSI  = getSequentialModelFromFI(modelFI);
modelASI = AdaptiveSequentialPressureTransportModel(modelSI.pressureModel, modelSI.transportModel, GF, G, rockF);

modelFIF  = TwoPhaseOilWaterModel(GF, rockF, fluid);
modelSIF  = getSequentialModelFromFI(modelFIF);

%%

time = 2*year;
rate = 2*sum(poreVolume(G, rock))/time;

xw = [0,0; l,l];

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

sW     = 0.0;
state0F = initResSol(GF, 100*barsa, [sW,1-sW]);

%%

modelASI.pressureModel.extraStateOutput = true;
[wsASI, stASI, rep] = simulateScheduleAD(state0, modelASI, schedule);

%%

[wsSI, stSI, rep] = simulateScheduleAD(state0, modelSI, schedule);

%%

[wsSIF, stSIF, rep] = simulateScheduleAD(state0F, modelSIF, scheduleF);

%%


close all
fig = figure('Position', [-2000, 0, 1500, 600]);
M = struct('cdata',[],'colormap',[]);
d = 10;

wcut = nan(numel(dt,1),3);
for sNo = 1:numel(stASI)
    wcut(sNo,:) = [wsSI{sNo}(2).wcut, wsASI{sNo}(2).wcut, wsSIF{sNo}(2).wcut];
end

clr = lines(3);
    
dt = cumsum(schedule.step.val)/day;

for sNo = 1:numel(stASI)
    
    subplot(1,3,1)
    plotCellData(stASI{sNo}.G, stASI{sNo}.s(:,1), 'edgec', 'none')
    plotGrid(stASI{sNo}.G, 'facec', 'none', 'edgealpha', 0.2);
    caxis([0,1]);
    axis equal tight
    
    subplot(1,3,2)
    plotCellData(G, stSI{sNo}.s(:,1), 'edgec', 'none')
    plotGrid(G, 'facec', 'none', 'edgealpha', 0.2);
    caxis([0,1]);
    axis equal tight
    
    colormap(jet)
    pause(0.1)
    
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

    
    rect = [d, d, fig.Position(3:4) - [d,d]];
%     rect = pos;
    M(sNo) = getframe(fig, rect);
    
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