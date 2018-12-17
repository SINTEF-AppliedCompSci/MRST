mrstModule add spe10 vem vemmech dg ad-core ad-props ad-blackoil ...
    blackoil-sequential gasinjection incomp msrsb

%%

[state0F, model, schedule]  = setupSPE10_AD('layers', 70);
GF    = model.G;
rockF = model.rock;
GF    = computeVEMGeometry(GF);
GF    = computeCellDimensions2(GF);

%%

close all
T = computeTrans(GF, rockF);
p = partitionMETIS(GF, T, 400);

GC = generateCoarseGrid(GF, p);
% plotGrid(GC)
plotGrid(GC, 10, 'facec', 'r')
axis equal tight

%%

% GC = generateCoarseGrid(GF, p);
GC = coarsenGeometry(GC);
GC = coarsenCellDimensions(GC);
GC = storeInteractionRegionCart(GC);

%%

fNo = 714;
clr = lines(3);
close all
figure; hold on
faces = G.faces.fconn(G.faces.connPos(fNo):G.faces.connPos(fNo+1)-1);
plotFaces(G.parent, faces, 'facealpha', 0.2);
set(gca, 'Zdir', 'normal')
axis equal tight; view(3)
xf = G.faces.centroids(fNo,:);
quiver3(xf(:,1), xf(:,2), xf(:,3), vec1(fNo,1), vec1(fNo,2), vec1(fNo,3), 20, 'linew', 2, 'color', clr(1,:));
quiver3(xf(:,1), xf(:,2), xf(:,3), vec2(fNo,1), vec2(fNo,2), vec2(fNo,3), 20, 'linew', 2, 'color', clr(1,:));

quiver3(xf(:,1), xf(:,2), xf(:,3), n(fNo, 1)   , n(fNo,2)   , n(fNo,3)   , 20 , 'linew', 2, 'color', clr(2,:));

%%

clr = lines(3);
close all
figure; hold on
fNo   = floor(rand*(GC.faces.num-1)) + 1;
faces = GC.faces.fconn(GC.faces.connPos(fNo):GC.faces.connPos(fNo+1)-1);
plotFaces(GF, faces, 'facealpha', 0.2);
axis equal tight

xf = GC.faces.centroids(fNo,:);
vec1 = GC.faces.coordSys{1}(fNo,:);
vec2 = GC.faces.coordSys{2}(fNo,:);
n    = GC.faces.normals(fNo,:)./GC.faces.areas(fNo);
quiver3(xf(:,1), xf(:,2), xf(:,3), vec1(1), vec1(2), vec1(3), 20, 'linew', 2, 'color', clr(1,:));
quiver3(xf(:,1), xf(:,2), xf(:,3), vec2(1), vec2(2), vec2(3), 20, 'linew', 2, 'color', clr(1,:));
quiver3(xf(:,1), xf(:,2), xf(:,3), n(1)   , n(2)   , n(3)   , 20 , 'linew', 2, 'color', clr(2,:));

xMax = GC.faces.xMax(fNo,:);
xMin = GC.faces.xMin(fNo,:);
plot3(xMax(:,1), xMax(:,2), xMax(:,3), '.');
plot3(xMin(:,1), xMin(:,2), xMin(:,3), '.');

x = GC.nodes.coords(GC.faces.nodes(GC.faces.nodePos(fNo):GC.faces.nodePos(fNo+1)-1),:);
xx = [x; x(1,:)];
patch(xx(:,1), xx(:,2), xx(:,3), 'r', 'facealpha', 0.2);
view(3)

%%

time = 2000*day;
bhp  = 50*barsa;
rate = 0.25*sum(poreVolume(GF, rockF))/time;

WF = [];
WF = verticalWell(WF, GF, rockF, 1, 1, []  , ...
                                   'type'  , 'rate', ...
                                   'val'   , rate  , ...
                                   'comp_i', [1,0] , ...
                                   'name'  , 'Inj' );
WF = verticalWell(WF, GF, rockF, 42, 220, [], ...
                                   'type'  , 'bhp' , ...
                                   'val'   , bhp   , ...
                                   'comp_i', [1,0] , ...
                                   'name'  , 'Prod');

%%

wc = GC.partition([WF.cells]);
mappings = getRefinementMappings(GC, GC, GF, wc);
G        = generateCoarseGrid(GF, mappings.newPartition);
G.cells.refined = mappings.refined;
G        = coarsenGeometry(G);
G        = coarsenCellDimensions(G);

%%

W = WF;
for wNo = 1:numel(W)
    W(wNo).cells = G.partition(W(wNo).cells);
end

schedule.control(1).W = W;
scheduleF = schedule;
scheduleF.control(1).W = WF;

%%

close all
plotToolbar(GF, rockF);
colormap(pink)
plotGrid(G, 'facec', 'none', 'edgec', [1,1,1]*0.25, 'linew', 1.5);
plotWell(G, W, 'color', 'm')
axis equal tight
ax = gca;
ax.Projection = 'perspective';

%%

modelFIF = TwoPhaseOilWaterModel(GF, rockF, model.fluid);
modelSIF = getSequentialModelFromFI(modelFIF);
modelASI = AdaptiveSequentialPressureTransportModel(modelSIF.pressureModel, modelSIF.transportModel, G);

modelFI = upscaleModelTPFA(modelFIF, G.partition);
% modelFI = upscaleHybridModelTPFA(modelFIF, G.partition);
modelSI = getSequentialModelFromFI(modelFI);

modelSIDG = getSequentialModelFromFI(modelFI);
GG = coarsenCellDimensions(modelSIDG.G);
modelSIDG.transportModel.G = GG;
[jt, ot, mt] = deal(Inf);
mt = 1e-3;
ot = 1e-3;
degree = 1;
discDG = DGDiscretization(modelSIDG.transportModel              , ...
                        'degree'               , degree       , ...
                        'basis'                , 'legendre'   , ...
                        'useUnstructCubature'  , false        , ...
                        'jumpTolerance'        , jt           , ...
                        'outTolerance'         , ot           , ...
                        'outLimiter'           , 'kill', ...
                        'meanTolerance'        , mt           , ...
                        'limitAfterConvergence', false        , ...
                        'plotLimiterProgress'  , false        );
transportModelDG = TransportOilWaterModelDG(G, modelFI.rock, fluid, ...
                                   'disc'    , discDG        , ...
                                   'dsMaxAbs', 0.2, ...
                                   'nonlinearTolerance', 1e-3);
modelSIDG.transportModel = transportModelDG;


modelSIFDG = getSequentialModelFromFI(modelFIF);
discDGF = DGDiscretization(modelSIFDG.transportModel          , ...
                        'degree'               , degree       , ...
                        'basis'                , 'legendre'   , ...
                        'useUnstructCubature'  , false        , ...
                        'jumpTolerance'        , jt           , ...
                        'outTolerance'         , ot           , ...
                        'outLimiter'           , 'kill', ...
                        'meanTolerance'        , mt           , ...
                        'limitAfterConvergence', false        , ...
                        'plotLimiterProgress'  , false        );
transportModelDGF = TransportOilWaterModelDG(GF, rockF, fluid, ...
                                   'disc'    , discDGF        , ...
                                   'dsMaxAbs', 0.2, ...
                                   'nonlinearTolerance', 1e-3);
modelSIFDG.transportModel = transportModelDGF;                        

%%

state0 = initResSol(G, state0F.pressure(1), state0F.s(1,:));
state0F.transportState        = state0;
state0F.transportState.degree = degree*ones(G.cells.num, 1);
state0F.transportState.G      = G;
state0F.transportState.pv     = modelASI.transportModel.operators.pv;

state0F.transportModel = modelASI.transportModel;
modelASI.plotProgress  = true;
close all
[wsASI, stASI, rep] = simulateScheduleAD(state0F, modelASI, scheduleF);

%%

state0F.transportModel = modelASI.transportModel;
[wsSIF, stSIF, rep] = simulateScheduleAD(state0F, modelSIF, scheduleF);

%%

[wsSI, stSI, rep] = simulateScheduleAD(state0, modelSI, schedule);

%%

% state0.degree = degree*ones(G.cells.num,1);
state0 = initResSol(G, state0F.pressure(1), state0F.s(1,:));
state0 = assignDofFromState(discDG, state0);
[wsSIDG, stSIDG, rep] = simulateScheduleAD(state0, modelSIDG, schedule);

%%

state0F = assignDofFromState(discDGF, state0F);
[wsSIFDG, stSIFDG, rep] = simulateScheduleAD(state0F, modelSIFDG, scheduleF);

%%

modelASI.storeGrids                     = true;
modelASI.plotProgress                   = true;
modelASI.pressureModel.extraStateOutput = true;

state0 = initResSol(G, state0F.pressure(1), state0F.s(1,:));
state0F.transportState        = state0;
state0F.transportState.degree = degree*ones(G.cells.num, 1);
state0F.transportState.G      = G;
state0F.transportState.pv     = modelASI.transportModel.operators.pv;
state0F.transportState        = discDG.assignDofFromState(state0F.transportState);
state0F.transportState        = discDG.updateDofPos(state0F.transportState);
state0F.transportModel        = modelASI.transportModel;

[wsASI, statesASI, rep] = simulateScheduleAD(state0F, modelASI, scheduleF);

%%

[ws, states, rep] = simulateScheduleAD(state0, model, schedule);

%%


%%

plotWellSols({wsSI, wsASI, wsSIF}, cumsum(scheduleF.step.val), 'datasetNames', {'Coarse', 'Adaptive', 'Fine'});

%%

close all;
pos = [-1000, 0, 600, 800];

%%

figure('position', pos)
plotToolbar(GF, stASI);
cmap = mrstColormap();
colormap(cmap);
axis equal tight; box on;

%%

figure('position', pos)
plotToolbar(GF, stSIF);
cmap = mrstColormap();
colormap(cmap);
axis equal tight; box on;

%%

figure('position', pos)
plotToolbar(G, stSI);
cmap = mrstColormap();
colormap(cmap);
axis equal tight; box on;

%%

dt = cumsum(schedule.step.val)/day;
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

for sNo = 1:numel(stASI)
    
    clf
    subplot(1,3,1)
    hold on
    plotCellData(modelASI.G, stASI{sNo}.s(:,1), 'edgec', 'none')
    plotGrid(stASI{sNo}.G, 'facecolor', 'none', 'edgealpha', 0.2);
    hold off
    caxis([0.2,0.8]);
    axis equal tight
    title('Adaptive')
    
    subplot(1,3,2)
    plotCellData(modelSIF.G, stSIF{sNo}.s(:,1), 'edgec', 'none')
    plotGrid(modelSIF.G, 'facecolor', 'none', 'edgealpha', 0.2);
    caxis([0.2,0.8]);
    axis equal tight
    title('Fine')
    
    colormap(cmap)
    
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
    legend({'Coarse', 'Adaptive', 'Fine'}, 'location', 'northwest');

    if 1
        rect = [d, d, fig.Position(3:4) - [d,d]];
        M(sNo) = getframe(fig, rect);
    end
    
    pause(0.1)
    
end

%%

pth = mrstPath('dg');
name = 'spe10-refinement-2';
duration = 10;
vo = VideoWriter(fullfile(pth, name));
vo.FrameRate = numel(stASI)/duration;
open(vo);

writeVideo(vo, M);

close(vo)

%%

close all

figure('position', [-2000, 0, 400, 800]);
plotToolbar(model.G, statesASI);

figure('position', [-2000, 0, 400, 800]);
plotToolbar(model.G, statesSIF);


%%

close all

fig = figure('Position', [-2000, 0, 1500, 1000]);

states = {statesDG, statesFV};
titles = {'dG(1)', 'FV'};

x = [250, 295];
for mNo = 1:numel(states)
    subplot(1,3,mNo+1)
    h(mNo) = plotCellData(G, states{mNo}{1}.s(:,1), 'edgec', 'none');
    colorbar('Location', 'southoutside');
    axis equal off
    text(x(mNo), 20, titles{mNo}, 'fontsize', 25, 'color', 'w'); 
    caxis([0,1])
end

subplot(1,3,1);
h(3) = plotCellData(G, states{1}{1}.degree, 'edgec', 'none');
colormap jet
caxis([0,1]);
colorbar('Location', 'southoutside');
axis equal off
text(150, 20, 'dG degree', 'fontsize', 25, 'color', 'w'); 

% set(fig, 'Units', 'pixels');
% pos = get(fig, 'Position');
% set(fig, 'Units', 'normalized');
M = struct('cdata',[],'colormap',[]);

for sNo = 1:numel(schedule.step.val)
    
    for mNo = 1:numel(states)
        h(mNo).CData = states{mNo}{sNo}.s(:,1);
    end
    
    h(3).CData = states{1}{sNo}.degree;
    
    pause(0.01);
    drawnow
%     dx = 10;
%     dy = 10;
    rect = [0, 0, fig.Position(3:4)];
%     rect = pos;
    M(sNo) = getframe(fig, rect);
    
end

%%

pth = mrstPath('dg');
name = 'spe10';
duration = 10;
vo = VideoWriter(fullfile(pth, name));
vo.FrameRate = numel(states{1})/duration;
open(vo);

writeVideo(vo, M);

close(vo)

